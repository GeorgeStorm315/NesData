/*
 * sn605p_adaptive_scurve_controller.c
 *
 * Reference implementation for the thesis:
 *   "Design and Validation of an Encoder-Synchronized Adaptive Jerk-Limited
 *    Speed-Control System for a Modified SN-605P Computerized Single-Cylinder
 *    Sock Knitting Machine"
 *
 * This is controller-independent research code, not the proprietary SN-605P
 * firmware. It is intended to run as a 1 kHz supervisory task above the
 * existing servo drive. In the preferred retrofit, omega_ref is written to the
 * drive through its documented speed-reference interface and the drive keeps
 * its original current/speed loops. The optional torque command is provided
 * for laboratory drives that expose a safe torque-mode interface.
 *
 * IMPORTANT SAFETY NOTE
 * ---------------------
 * This software is not a safety function. Emergency stop, guard interlock,
 * safe torque off, overtravel protection, and drive protection must remain
 * hardwired or implemented in a safety-rated system according to the machine
 * risk assessment and applicable standards. Commission first with the needle
 * cylinder unloaded and with conservative limits.
 *
 * Build demonstration:
 *   cc -std=c11 -O2 -Wall -Wextra -pedantic \
 *      -DSN605P_CONTROLLER_DEMO sn605p_adaptive_scurve_controller.c \
 *      -lm -o sn605p_demo
 *
 * The demo writes comma-separated simulation data to standard output.
 *
 * SPDX-License-Identifier: MIT
 */

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SN_ARRAY_LEN(a) (sizeof(a) / sizeof((a)[0]))

typedef enum {
    SN_SECTION_IDLE = 0,
    SN_SECTION_CUFF,
    SN_SECTION_LEG,
    SN_SECTION_ANKLE,
    SN_SECTION_HEEL_ENTRY,
    SN_SECTION_HEEL_RECIPROCATION,
    SN_SECTION_FOOT,
    SN_SECTION_TOE_ENTRY,
    SN_SECTION_TOE_RECIPROCATION,
    SN_SECTION_RUNOUT,
    SN_SECTION_FAULT,
    SN_SECTION_COUNT
} sn_section_t;

typedef enum {
    SN_MODE_DISABLED = 0,
    SN_MODE_READY,
    SN_MODE_RUNNING,
    SN_MODE_CONTROLLED_STOP,
    SN_MODE_FAULT_LATCHED
} sn_mode_t;

typedef struct {
    /* Supervisory sample time. */
    double dt_s;

    /* Identified equivalent mechanical parameters at the cylinder shaft. */
    double inertia_kgm2;
    double viscous_friction_nms;

    /* Reference-generator limits. */
    double accel_limit_rad_s2;
    double decel_limit_rad_s2;
    double jerk_limit_rad_s3;
    double reference_error_gain_s_inv;

    /* Load-observer settings. */
    double observer_bandwidth_rad_s;
    double accel_filter_bandwidth_rad_s;
    double rated_torque_nm;

    /* Adaptive-limit floors and utilization threshold. */
    double minimum_limit_scale;
    double utilization_soft_limit;

    /* Optional PI torque loop for a laboratory torque-mode drive. */
    double speed_kp_nm_per_rad_s;
    double speed_ki_nm_per_rad;
    double torque_limit_nm;
    double antiwindup_gain_s_inv;

    /* Once-per-revolution index correction and actuator delay. */
    double index_correction_gain;
    double max_index_correction_rad;
    double actuator_delay_s;
    double calibrated_phase_offset_rad;

    /* Monitoring thresholds. */
    double overspeed_rad_s;
    double jam_speed_error_rad_s;
    double jam_load_threshold_nm;
    double jam_persistence_s;
} sn_controller_config_t;

typedef struct {
    /* Pattern engine supplies the required signed speed and current section. */
    sn_section_t section;
    double requested_speed_rad_s;

    /* Encoder and drive feedback. */
    double encoder_angle_rad;
    double measured_speed_rad_s;
    double applied_torque_nm;
    double torque_or_current_utilization; /* normalized 0.0 ... 1.0+ */
    bool index_pulse;

    /* Machine permissives and stops. */
    bool enable_request;
    bool drive_ready;
    bool yarn_break;
    bool guard_open;
    bool emergency_stop;
    bool reset_fault;
} sn_controller_input_t;

typedef struct {
    /* Command to a standard speed-mode servo retrofit. */
    double speed_reference_rad_s;
    double accel_reference_rad_s2;

    /* Optional command for a torque-mode research drive. */
    double torque_reference_nm;

    /* Diagnostics for logging and commissioning. */
    double estimated_load_torque_nm;
    double adaptive_limit_scale;
    double phase_correction_rad;
    sn_mode_t mode;
    uint32_t fault_bits;
} sn_controller_output_t;

typedef struct {
    sn_mode_t mode;
    double omega_ref_rad_s;
    double alpha_ref_rad_s2;
    double previous_measured_speed_rad_s;
    double filtered_accel_rad_s2;
    double estimated_load_torque_nm;
    double speed_integral_rad;
    double phase_correction_rad;
    double jam_timer_s;
    uint32_t fault_bits;
    bool initialized;
} sn_controller_state_t;

enum {
    SN_FAULT_ESTOP       = 1u << 0,
    SN_FAULT_GUARD       = 1u << 1,
    SN_FAULT_DRIVE       = 1u << 2,
    SN_FAULT_YARN        = 1u << 3,
    SN_FAULT_OVERSPEED   = 1u << 4,
    SN_FAULT_JAM         = 1u << 5,
    SN_FAULT_BAD_SECTION = 1u << 6,
    SN_FAULT_BAD_CONFIG  = 1u << 7
};

typedef struct {
    double speed_cap_rad_s;
    double accel_scale;
    double jerk_scale;
} sn_section_limits_t;

static const sn_section_limits_t k_section_limits[SN_SECTION_COUNT] = {
    [SN_SECTION_IDLE]                = {0.0,  0.50, 0.40},
    [SN_SECTION_CUFF]                = {24.1, 0.65, 0.55}, /* 230 rpm */
    [SN_SECTION_LEG]                 = {33.5, 1.00, 1.00}, /* 320 rpm */
    [SN_SECTION_ANKLE]               = {27.2, 0.80, 0.70}, /* 260 rpm */
    [SN_SECTION_HEEL_ENTRY]          = {18.8, 0.55, 0.45}, /* 180 rpm */
    [SN_SECTION_HEEL_RECIPROCATION]  = {13.6, 0.45, 0.35}, /* 130 rpm */
    [SN_SECTION_FOOT]                = {29.3, 0.90, 0.85}, /* 280 rpm */
    [SN_SECTION_TOE_ENTRY]           = {17.8, 0.55, 0.45}, /* 170 rpm */
    [SN_SECTION_TOE_RECIPROCATION]   = {12.6, 0.42, 0.32}, /* 120 rpm */
    [SN_SECTION_RUNOUT]              = {10.5, 0.45, 0.35}, /* 100 rpm */
    [SN_SECTION_FAULT]               = {0.0,  0.35, 0.25}
};

static double sn_clamp(double x, double lower, double upper)
{
    return (x < lower) ? lower : ((x > upper) ? upper : x);
}

static double sn_abs(double x)
{
    return (x < 0.0) ? -x : x;
}

static double sn_wrap_0_2pi(double angle_rad)
{
    double wrapped = fmod(angle_rad, 2.0 * M_PI);
    if (wrapped < 0.0) {
        wrapped += 2.0 * M_PI;
    }
    return wrapped;
}

static double sn_lowpass_alpha(double bandwidth_rad_s, double dt_s)
{
    const double x = bandwidth_rad_s * dt_s;
    return x / (1.0 + x);
}

static bool sn_config_is_valid(const sn_controller_config_t *cfg)
{
    return cfg != NULL &&
           cfg->dt_s > 0.0 &&
           cfg->inertia_kgm2 > 0.0 &&
           cfg->accel_limit_rad_s2 > 0.0 &&
           cfg->decel_limit_rad_s2 > 0.0 &&
           cfg->jerk_limit_rad_s3 > 0.0 &&
           cfg->rated_torque_nm > 0.0 &&
           cfg->torque_limit_nm > 0.0 &&
           cfg->minimum_limit_scale > 0.0 &&
           cfg->minimum_limit_scale <= 1.0;
}

sn_controller_config_t sn_controller_default_config(void)
{
    sn_controller_config_t cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.dt_s                         = 0.001;
    cfg.inertia_kgm2                 = 0.018;
    cfg.viscous_friction_nms         = 0.030;
    cfg.accel_limit_rad_s2           = 85.0;
    cfg.decel_limit_rad_s2           = 105.0;
    cfg.jerk_limit_rad_s3            = 950.0;
    cfg.reference_error_gain_s_inv   = 7.0;
    cfg.observer_bandwidth_rad_s      = 18.0;
    cfg.accel_filter_bandwidth_rad_s  = 55.0;
    cfg.rated_torque_nm               = 7.0;
    cfg.minimum_limit_scale           = 0.35;
    cfg.utilization_soft_limit        = 0.72;
    cfg.speed_kp_nm_per_rad_s         = 0.65;
    cfg.speed_ki_nm_per_rad           = 8.0;
    cfg.torque_limit_nm               = 8.0;
    cfg.antiwindup_gain_s_inv         = 12.0;
    cfg.index_correction_gain         = 0.08;
    cfg.max_index_correction_rad      = 0.012;
    cfg.actuator_delay_s              = 0.0045;
    cfg.calibrated_phase_offset_rad   = 0.0;
    cfg.overspeed_rad_s               = 38.0;
    cfg.jam_speed_error_rad_s         = 8.0;
    cfg.jam_load_threshold_nm         = 5.5;
    cfg.jam_persistence_s             = 0.120;
    return cfg;
}

void sn_controller_init(sn_controller_state_t *state)
{
    if (state != NULL) {
        memset(state, 0, sizeof(*state));
        state->mode = SN_MODE_DISABLED;
        state->initialized = true;
    }
}

/*
 * Compute a firing angle that anticipates a known actuator/pneumatic delay.
 * theta_sync is a slowly learned once-per-revolution correction; it must not
 * be changed discontinuously while knitting.
 */
double sn_compensated_firing_angle(double target_angle_rad,
                                   double omega_rad_s,
                                   double alpha_rad_s2,
                                   double actuator_delay_s,
                                   double calibrated_offset_rad,
                                   double theta_sync_rad)
{
    const double predicted_motion =
        omega_rad_s * actuator_delay_s +
        0.5 * alpha_rad_s2 * actuator_delay_s * actuator_delay_s;
    return sn_wrap_0_2pi(target_angle_rad - predicted_motion +
                         calibrated_offset_rad + theta_sync_rad);
}

static uint32_t sn_monitor_faults(const sn_controller_config_t *cfg,
                                  sn_controller_state_t *state,
                                  const sn_controller_input_t *in)
{
    uint32_t faults = state->fault_bits;

    if (in->emergency_stop) {
        faults |= SN_FAULT_ESTOP;
    }
    if (in->guard_open) {
        faults |= SN_FAULT_GUARD;
    }
    if (!in->drive_ready) {
        faults |= SN_FAULT_DRIVE;
    }
    if (in->yarn_break) {
        faults |= SN_FAULT_YARN;
    }
    if (sn_abs(in->measured_speed_rad_s) > cfg->overspeed_rad_s) {
        faults |= SN_FAULT_OVERSPEED;
    }
    if (in->section < SN_SECTION_IDLE || in->section >= SN_SECTION_COUNT) {
        faults |= SN_FAULT_BAD_SECTION;
    }

    const bool jam_condition =
        sn_abs(state->omega_ref_rad_s - in->measured_speed_rad_s) >
            cfg->jam_speed_error_rad_s &&
        sn_abs(state->estimated_load_torque_nm) >
            cfg->jam_load_threshold_nm;
    if (jam_condition) {
        state->jam_timer_s += cfg->dt_s;
    } else {
        state->jam_timer_s = fmax(0.0, state->jam_timer_s - 2.0 * cfg->dt_s);
    }
    if (state->jam_timer_s >= cfg->jam_persistence_s) {
        faults |= SN_FAULT_JAM;
    }

    if (in->reset_fault && !in->emergency_stop && !in->guard_open &&
        in->drive_ready && !in->yarn_break &&
        sn_abs(in->measured_speed_rad_s) < 0.5) {
        faults = 0u;
        state->jam_timer_s = 0.0;
    }
    return faults;
}

static double sn_adaptive_scale(const sn_controller_config_t *cfg,
                                const sn_controller_input_t *in,
                                double estimated_load_nm)
{
    const double load_ratio = sn_clamp(sn_abs(estimated_load_nm) /
                                       cfg->rated_torque_nm, 0.0, 1.5);
    const double utilization = sn_clamp(in->torque_or_current_utilization,
                                        0.0, 1.5);
    const double load_penalty = 0.55 * load_ratio;
    const double current_penalty =
        (utilization > cfg->utilization_soft_limit)
            ? 1.5 * (utilization - cfg->utilization_soft_limit)
            : 0.0;
    return sn_clamp(1.0 - load_penalty - current_penalty,
                    cfg->minimum_limit_scale, 1.0);
}

static void sn_update_observer(const sn_controller_config_t *cfg,
                               sn_controller_state_t *state,
                               const sn_controller_input_t *in)
{
    const double raw_accel =
        (in->measured_speed_rad_s - state->previous_measured_speed_rad_s) /
        cfg->dt_s;
    const double accel_alpha = sn_lowpass_alpha(
        cfg->accel_filter_bandwidth_rad_s, cfg->dt_s);
    state->filtered_accel_rad_s2 += accel_alpha *
        (raw_accel - state->filtered_accel_rad_s2);

    const double raw_load_nm =
        in->applied_torque_nm -
        cfg->inertia_kgm2 * state->filtered_accel_rad_s2 -
        cfg->viscous_friction_nms * in->measured_speed_rad_s;
    const double observer_alpha = sn_lowpass_alpha(
        cfg->observer_bandwidth_rad_s, cfg->dt_s);
    state->estimated_load_torque_nm += observer_alpha *
        (raw_load_nm - state->estimated_load_torque_nm);
    state->previous_measured_speed_rad_s = in->measured_speed_rad_s;
}

static void sn_update_index_correction(const sn_controller_config_t *cfg,
                                       sn_controller_state_t *state,
                                       const sn_controller_input_t *in)
{
    if (!in->index_pulse) {
        return;
    }
    /* At the physical index, wrapped encoder angle should be zero. */
    double phase_error = sn_wrap_0_2pi(in->encoder_angle_rad);
    if (phase_error > M_PI) {
        phase_error -= 2.0 * M_PI;
    }
    const double correction = sn_clamp(
        -cfg->index_correction_gain * phase_error,
        -cfg->max_index_correction_rad,
         cfg->max_index_correction_rad);
    state->phase_correction_rad = sn_clamp(
        state->phase_correction_rad + correction,
        -cfg->max_index_correction_rad,
         cfg->max_index_correction_rad);
}

static void sn_update_reference(const sn_controller_config_t *cfg,
                                sn_controller_state_t *state,
                                const sn_controller_input_t *in,
                                double adaptive_scale)
{
    sn_section_t section = in->section;
    if (section < SN_SECTION_IDLE || section >= SN_SECTION_COUNT) {
        section = SN_SECTION_FAULT;
    }
    const sn_section_limits_t lim = k_section_limits[section];

    double target = in->requested_speed_rad_s;
    target = sn_clamp(target, -lim.speed_cap_rad_s, lim.speed_cap_rad_s);
    if (state->mode != SN_MODE_RUNNING) {
        target = 0.0;
    }

    const double section_scale = lim.accel_scale * adaptive_scale;
    const double accel_limit = cfg->accel_limit_rad_s2 * section_scale;
    const double decel_limit = cfg->decel_limit_rad_s2 * section_scale;
    const double jerk_limit = cfg->jerk_limit_rad_s3 *
                              lim.jerk_scale * adaptive_scale;
    const double speed_error = target - state->omega_ref_rad_s;

    /* Online jerk-limited reference generator. The proportional mapping from
     * speed error to desired acceleration yields the seven-segment S-shaped
     * response when limits remain constant, but it can also retarget safely
     * between consecutive courses. */
    double desired_accel = cfg->reference_error_gain_s_inv * speed_error;
    desired_accel = sn_clamp(desired_accel, -decel_limit, accel_limit);

    const double max_accel_step = jerk_limit * cfg->dt_s;
    const double accel_delta = sn_clamp(
        desired_accel - state->alpha_ref_rad_s2,
        -max_accel_step, max_accel_step);
    state->alpha_ref_rad_s2 += accel_delta;

    double next_speed = state->omega_ref_rad_s +
                        state->alpha_ref_rad_s2 * cfg->dt_s;

    /* Prevent a discrete integration step from crossing a stationary target. */
    if ((speed_error > 0.0 && next_speed > target) ||
        (speed_error < 0.0 && next_speed < target)) {
        next_speed = target;
        state->alpha_ref_rad_s2 = 0.0;
    }
    if (sn_abs(speed_error) < 1.0e-6 &&
        sn_abs(state->alpha_ref_rad_s2) < max_accel_step) {
        state->alpha_ref_rad_s2 = 0.0;
        next_speed = target;
    }
    state->omega_ref_rad_s = next_speed;
}

static double sn_optional_torque_command(const sn_controller_config_t *cfg,
                                         sn_controller_state_t *state,
                                         const sn_controller_input_t *in)
{
    const double error = state->omega_ref_rad_s - in->measured_speed_rad_s;
    const double feedforward =
        cfg->inertia_kgm2 * state->alpha_ref_rad_s2 +
        cfg->viscous_friction_nms * state->omega_ref_rad_s +
        state->estimated_load_torque_nm;
    const double unsaturated =
        feedforward + cfg->speed_kp_nm_per_rad_s * error +
        cfg->speed_ki_nm_per_rad * state->speed_integral_rad;
    const double saturated = sn_clamp(unsaturated,
                                      -cfg->torque_limit_nm,
                                       cfg->torque_limit_nm);

    /* Back-calculation anti-windup. */
    const double aw = cfg->antiwindup_gain_s_inv *
                      (saturated - unsaturated) /
                      fmax(cfg->speed_ki_nm_per_rad, 1.0e-9);
    state->speed_integral_rad += (error + aw) * cfg->dt_s;
    return saturated;
}

sn_controller_output_t sn_controller_step(const sn_controller_config_t *cfg,
                                          sn_controller_state_t *state,
                                          const sn_controller_input_t *in)
{
    sn_controller_output_t out;
    memset(&out, 0, sizeof(out));

    if (cfg == NULL || state == NULL || in == NULL ||
        !state->initialized || !sn_config_is_valid(cfg)) {
        out.mode = SN_MODE_FAULT_LATCHED;
        out.fault_bits = SN_FAULT_BAD_CONFIG;
        return out;
    }

    sn_update_observer(cfg, state, in);
    sn_update_index_correction(cfg, state, in);
    state->fault_bits = sn_monitor_faults(cfg, state, in);

    if (state->fault_bits != 0u) {
        state->mode = SN_MODE_FAULT_LATCHED;
    } else if (!in->enable_request) {
        state->mode = (sn_abs(in->measured_speed_rad_s) > 0.5)
                        ? SN_MODE_CONTROLLED_STOP : SN_MODE_READY;
    } else if (in->drive_ready) {
        state->mode = SN_MODE_RUNNING;
    } else {
        state->mode = SN_MODE_DISABLED;
    }

    const double adaptive_scale = sn_adaptive_scale(
        cfg, in, state->estimated_load_torque_nm);
    sn_update_reference(cfg, state, in, adaptive_scale);

    out.speed_reference_rad_s = state->omega_ref_rad_s;
    out.accel_reference_rad_s2 = state->alpha_ref_rad_s2;
    out.torque_reference_nm = sn_optional_torque_command(cfg, state, in);
    out.estimated_load_torque_nm = state->estimated_load_torque_nm;
    out.adaptive_limit_scale = adaptive_scale;
    out.phase_correction_rad = state->phase_correction_rad;
    out.mode = state->mode;
    out.fault_bits = state->fault_bits;
    return out;
}

const char *sn_section_name(sn_section_t section)
{
    static const char *const names[SN_SECTION_COUNT] = {
        "idle", "cuff", "leg", "ankle", "heel_entry",
        "heel_recip", "foot", "toe_entry", "toe_recip",
        "runout", "fault"
    };
    return (section >= SN_SECTION_IDLE && section < SN_SECTION_COUNT)
        ? names[section] : "invalid";
}

#ifdef SN605P_CONTROLLER_DEMO

typedef struct {
    sn_section_t section;
    double command_rad_s;
} demo_pattern_t;

static demo_pattern_t demo_pattern(double t_s)
{
    if (t_s < 0.5) {
        return (demo_pattern_t){SN_SECTION_IDLE, 0.0};
    }
    if (t_s < 3.0) {
        return (demo_pattern_t){SN_SECTION_LEG, 300.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 3.7) {
        return (demo_pattern_t){SN_SECTION_HEEL_ENTRY,
                                150.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 7.3) {
        const int half_cycle = (int)((t_s - 3.7) / 0.45);
        const double direction = (half_cycle % 2 == 0) ? 1.0 : -1.0;
        return (demo_pattern_t){SN_SECTION_HEEL_RECIPROCATION,
                                direction * 115.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 10.0) {
        return (demo_pattern_t){SN_SECTION_FOOT,
                                270.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 10.7) {
        return (demo_pattern_t){SN_SECTION_TOE_ENTRY,
                                145.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 13.4) {
        const int half_cycle = (int)((t_s - 10.7) / 0.40);
        const double direction = (half_cycle % 2 == 0) ? 1.0 : -1.0;
        return (demo_pattern_t){SN_SECTION_TOE_RECIPROCATION,
                                direction * 105.0 * 2.0 * M_PI / 60.0};
    }
    if (t_s < 14.5) {
        return (demo_pattern_t){SN_SECTION_RUNOUT,
                                80.0 * 2.0 * M_PI / 60.0};
    }
    return (demo_pattern_t){SN_SECTION_IDLE, 0.0};
}

int main(void)
{
    sn_controller_config_t cfg = sn_controller_default_config();
    sn_controller_state_t controller;
    sn_controller_init(&controller);

    double plant_speed = 0.0;
    double plant_angle = 0.0;
    double applied_torque = 0.0;
    const double plant_inertia = 0.020;
    const double plant_friction = 0.034;
    const double duration_s = 16.0;

    puts("time_s,section,command_rpm,reference_rpm,measured_rpm,"
         "accel_ref_rad_s2,torque_nm,load_est_nm,limit_scale,phase_fire_deg");

    for (size_t k = 0; k < (size_t)(duration_s / cfg.dt_s); ++k) {
        const double t = (double)k * cfg.dt_s;
        const demo_pattern_t pattern = demo_pattern(t);
        const double actual_load = 1.1 + 0.20 * sin(3.0 * plant_angle) +
            ((pattern.section == SN_SECTION_HEEL_RECIPROCATION ||
              pattern.section == SN_SECTION_TOE_RECIPROCATION) ? 1.35 : 0.0) +
            ((t > 8.2 && t < 8.8) ? 2.2 : 0.0);

        const double previous_angle = plant_angle;
        plant_angle = sn_wrap_0_2pi(plant_angle + plant_speed * cfg.dt_s);
        const bool index_pulse = plant_angle < previous_angle && plant_speed >= 0.0;

        sn_controller_input_t in;
        memset(&in, 0, sizeof(in));
        in.section = pattern.section;
        in.requested_speed_rad_s = pattern.command_rad_s;
        in.encoder_angle_rad = plant_angle;
        in.measured_speed_rad_s = plant_speed;
        in.applied_torque_nm = applied_torque;
        in.torque_or_current_utilization = sn_abs(applied_torque) /
                                           cfg.torque_limit_nm;
        in.index_pulse = index_pulse;
        in.enable_request = t >= 0.5 && t < 15.0;
        in.drive_ready = true;

        const sn_controller_output_t out = sn_controller_step(&cfg, &controller, &in);
        applied_torque = out.torque_reference_nm;
        const double load_direction =
            (sn_abs(plant_speed) > 0.02) ? copysign(1.0, plant_speed) :
            ((sn_abs(applied_torque) > actual_load) ?
                copysign(1.0, applied_torque) : 0.0);
        const double acceleration = (applied_torque -
            plant_friction * plant_speed - actual_load * load_direction) /
            plant_inertia;
        plant_speed += acceleration * cfg.dt_s;

        const double firing = sn_compensated_firing_angle(
            210.0 * M_PI / 180.0, out.speed_reference_rad_s,
            out.accel_reference_rad_s2, cfg.actuator_delay_s,
            cfg.calibrated_phase_offset_rad, out.phase_correction_rad);

        if (k % 10u == 0u) {
            printf("%.3f,%s,%.3f,%.3f,%.3f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                   t, sn_section_name(pattern.section),
                   pattern.command_rad_s * 60.0 / (2.0 * M_PI),
                   out.speed_reference_rad_s * 60.0 / (2.0 * M_PI),
                   plant_speed * 60.0 / (2.0 * M_PI),
                   out.accel_reference_rad_s2, applied_torque,
                   out.estimated_load_torque_nm, out.adaptive_limit_scale,
                   firing * 180.0 / M_PI);
        }
    }
    return 0;
}

#endif /* SN605P_CONTROLLER_DEMO */
