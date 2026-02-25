/**
 * @file advection.cpp
 * @brief Implementation for the advection module.
 *
 * Provides executable logic for the advection runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/advection subsystem.
 */

#include "advection.hpp"
#include "advection_base.hpp"
#include "simulation.hpp"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <string>



namespace
{

/**
 * @brief Struct for the advection performance totals.
 */
struct AdvectionPerfTotals
{
    uint64_t scalar_calls = 0;
    uint64_t thermo_calls = 0;
    uint64_t tracer_calls = 0;
    double scalar_total_s = 0.0;
    double thermo_total_s = 0.0;
    double tracer_total_s = 0.0;
    double r_s = 0.0;
    double theta_s = 0.0;
    double z_s = 0.0;
    double diffusion_s = 0.0;
};

/**
 * @brief Aggregated advection timing counters for runtime diagnostics.
 */
AdvectionPerfTotals g_advection_perf_totals;

/**
 * @brief Returns whether advection performance timing is enabled.
 */
inline bool perf_enabled() {return global_perf_timing_enabled;}


/**
 * @brief Resizes a field to the active domain dimensions when needed.
 */
inline void ensure_field_shape(Field3D& field)
{

    if (field.size_r() != NR || field.size_th() != NTH || field.size_z() != NZ)
    {
        field.resize(NR, NTH, NZ, 0.0f);
    }
}

/**
 * @brief Copies one 3D field into another with shape validation.
 */
inline void copy_field(const Field3D& src, Field3D& dst)
{
    ensure_field_shape(dst);

    if (src.size() == 0){return;}

    std::memcpy(dst.data(), src.data(), src.size() * sizeof(float));
}

/**
 * @brief Computes the row-major flat index for a 3D cell.
 */
inline size_t idx3(const int i, const int j, const int k)
{
    return (static_cast<size_t>(i) * static_cast<size_t>(NTH) + static_cast<size_t>(j)) *
           static_cast<size_t>(NZ) + static_cast<size_t>(k);
}

/**
 * @brief Copies radial and vertical boundary values for cylindrical domains.
 */
inline void copy_cylindrical_boundaries(const Field3D& src, Field3D& dst)
{
    if (src.size() == 0){return;}

    const float* src_data = src.data();
    float* dst_data = dst.data();

    const size_t plane_size = static_cast<size_t>(NTH) * static_cast<size_t>(NZ);

    std::memcpy(dst_data, src_data, plane_size * sizeof(float));

    if (NR > 1)
    {
        const size_t last_plane = (static_cast<size_t>(NR) - 1) * plane_size;
        std::memcpy(dst_data + last_plane, src_data + last_plane, plane_size * sizeof(float));
    }

    if (NZ > 0 && NR > 2)
    {

        #pragma omp parallel for collapse(2) 

        for (int i = 1; i < NR - 1; ++i)
        {
            for (int j = 0; j < NTH; ++j)
            {
                const size_t k0 = idx3(i, j, 0);
                dst_data[k0] = src_data[k0];

                if (NZ > 1)
                {
                    const size_t ktop = idx3(i, j, NZ - 1);
                    dst_data[ktop] = src_data[ktop];
                }
            }
        }
    }
}

/**
 * @brief Returns true when any field entry exceeds a small threshold.
 */
inline bool field_has_signal(const Field3D& field, float eps = 1.0e-12f)
{
    const float* data = field.data(); 
    const size_t count = field.size();

    for (size_t idx = 0; idx < count; ++idx)
    {
        if (std::fabs(data[idx]) > eps)
        {
            return true;
        }
    }
    return false;
}

/**
 * @brief Returns a lowercase copy of an input string.
 */
std::string lower_copy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

/**
 * @brief Checks whether the runtime should use numerics vertical transport.
 */
bool use_numerics_vertical_advection()
{
    if (!advection_scheme){return false;}

    const std::string scheme_id = lower_copy(advection_scheme->name());
    return scheme_id == "tvd" || scheme_id == "weno5";
}

/**
 * @brief Logs the selected runtime advection path once per execution.
 */
void log_runtime_advection_path_once()
{
    static bool logged = false;
    if (logged || !log_normal_enabled()){return;}
    logged = true;

    if (use_numerics_vertical_advection())
    {
        std::cout << "[ADVECTION] Runtime path: src/advection directional split (r/theta) + "
                  << "src/numerics/advection/" << global_advection_config.scheme_id
                  << " (vertical tendency)" << std::endl;
    }
    else
    {
        std::cout << "[ADVECTION] Runtime path: src/advection directional split kernels "
                  << "(numerics.advection schemes are initialized but not selected for vertical transport)." << std::endl;
    }
}

}

/**
 * @brief Applies explicit scalar diffusion on the interior domain.
 */
static void apply_diffusion_kernel(const Field3D& src, Field3D& dst, double dt, double kappa)
{
    ensure_field_shape(dst);

    if (kappa <= 0.0)
    {
        copy_field(src, dst);
        return;
    }

    const float* src_data = src.data();
    float* dst_data = dst.data();

    copy_cylindrical_boundaries(src, dst);
    const double inv_dr2 = 1.0 / (dr * dr);
    const double inv_dz2 = 1.0 / (dz * dz);

    #pragma omp parallel for collapse(2) 

    for (int i = 1; i < NR - 1; ++i)
    {
        const double r = i * dr + 1.0e-6;

        const double inv_r2_dtheta2 = 1.0 / (dtheta * dtheta * r * r);

        for (int j = 0; j < NTH; ++j)
        {
            const int j_prev = (j - 1 + NTH) % NTH;
            const int j_next = (j + 1) % NTH;

            for (int k = 1; k < NZ - 1; ++k)
            {
                const size_t c = idx3(i, j, k);
                const size_t ip = idx3(i + 1, j, k);
                const size_t im = idx3(i - 1, j, k);
                const size_t kp = idx3(i, j, k + 1);
                const size_t km = idx3(i, j, k - 1);
                const size_t jp = idx3(i, j_next, k);
                const size_t jm = idx3(i, j_prev, k);

                const double q = static_cast<double>(src_data[c]);
                const double lap_s =
                    (static_cast<double>(src_data[ip]) - 2.0 * q + static_cast<double>(src_data[im])) * inv_dr2 +
                    (static_cast<double>(src_data[kp]) - 2.0 * q + static_cast<double>(src_data[km])) * inv_dz2 +
                    (static_cast<double>(src_data[jp]) - 2.0 * q + static_cast<double>(src_data[jm])) * inv_r2_dtheta2;

                dst_data[c] = static_cast<float>(q + dt * kappa * lap_s);
            }
        }
    }
}

/**
 * @brief Advects a scalar in the radial direction with first-order upwinding.
 */
static void advect_scalar_1d_r_kernel(const Field3D& src, Field3D& dst, double dt, double)
{
    ensure_field_shape(dst);

    copy_cylindrical_boundaries(src, dst);

    const float* src_data = src.data();
    const float* u_data = u.data();
    float* dst_data = dst.data();

    
    #pragma omp parallel for collapse(2)

    for (int i = 1; i < NR - 1; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            for (int k = 1; k < NZ - 1; ++k)
            {
                const size_t c = idx3(i, j, k);
                const size_t im = idx3(i - 1, j, k);
                const size_t ip = idx3(i + 1, j, k);

                const double u_val = static_cast<double>(u_data[c]);
                double dq_dr = 0.0;

                if (u_val > 0.0)
                {
                    dq_dr = (static_cast<double>(src_data[c]) - static_cast<double>(src_data[im])) / dr;
                }
                else if (u_val < 0.0)
                {
                    dq_dr = (static_cast<double>(src_data[ip]) - static_cast<double>(src_data[c])) / dr;
                }
                const double q = static_cast<double>(src_data[c]);
                dst_data[c] = static_cast<float>(q - dt * u_val * dq_dr);
            }
        }
    }
}

/**
 * @brief Advects a scalar in the azimuthal direction with first-order upwinding.
 */
static void advect_scalar_1d_theta_kernel(const Field3D& src, Field3D& dst, double dt, double)
{
    ensure_field_shape(dst);

    copy_cylindrical_boundaries(src, dst);

    const float* src_data = src.data();
    const float* v_data = v_theta.data();
    float* dst_data = dst.data();

    
    #pragma omp parallel for collapse(2)

    for (int i = 1; i < NR - 1; ++i)
    {
        const double r = i * dr + 1.0e-6;

        for (int j = 0; j < NTH; ++j)
        {
            const int j_prev = (j - 1 + NTH) % NTH;
            const int j_next = (j + 1) % NTH;

            for (int k = 1; k < NZ - 1; ++k)
            {
                const size_t c = idx3(i, j, k);
                const size_t jm = idx3(i, j_prev, k);
                const size_t jp = idx3(i, j_next, k);

                const double v_val = static_cast<double>(v_data[c]);
                double dq_dtheta = 0.0;

                if (v_val > 0.0)
                {
                    dq_dtheta = (static_cast<double>(src_data[c]) - static_cast<double>(src_data[jm])) / dtheta;
                }
                else if (v_val < 0.0)
                {
                    dq_dtheta = (static_cast<double>(src_data[jp]) - static_cast<double>(src_data[c])) / dtheta;
                }

                const double dqdt = -(v_val / r) * dq_dtheta;
                dst_data[c] = static_cast<float>(static_cast<double>(src_data[c]) + dt * dqdt);
            }
        }
    }
}


/**
 * @brief Advects a scalar in the vertical direction with first-order upwinding.
 */
static void advect_scalar_1d_z_kernel(const Field3D& src, Field3D& dst, double dt, double)
{
    ensure_field_shape(dst);
    copy_cylindrical_boundaries(src, dst);

    const float* src_data = src.data();
    const float* w_data = w.data();
    float* dst_data = dst.data();

    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NR - 1; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            for (int k = 1; k < NZ - 1; ++k)
            {
                const size_t c = idx3(i, j, k);
                const size_t km = idx3(i, j, k - 1);
                const size_t kp = idx3(i, j, k + 1);

                const double w_val = static_cast<double>(w_data[c]);
                double dq_dz = 0.0;
                if (w_val > 0.0)
                {
                    dq_dz = (static_cast<double>(src_data[c]) - static_cast<double>(src_data[km])) / dz;
                }
                else if (w_val < 0.0)
                {
                    dq_dz = (static_cast<double>(src_data[kp]) - static_cast<double>(src_data[c])) / dz;
                }
                const double q = static_cast<double>(src_data[c]);
                dst_data[c] = static_cast<float>(q - dt * w_val * dq_dz);
            }
        }
    }
}

/**
 * @brief Applies numerics-module vertical advection and falls back on failure.
 */
static void advect_scalar_1d_z_numerics_kernel(const Field3D& src, Field3D& dst, double dt, double)
{
    ensure_field_shape(dst);
    copy_field(src, dst);
    copy_cylindrical_boundaries(src, dst);

    if (!advection_scheme || !use_numerics_vertical_advection())
    {
        return;
    }

    static AdvectionTendencies tendencies;
    AdvectionStateView state{};
    state.u = &u;
    state.v = &v_theta;
    state.w = &w;
    state.q = &src;
    state.rho = &rho;
    state.grid = &global_grid_metrics;
    AdvectionConfig runtime_cfg = global_advection_config;
    runtime_cfg.positivity_dt = std::max(std::abs(dt), 1.0e-12);

    try
    {
        advection_scheme->compute_flux_divergence(runtime_cfg, state, tendencies, nullptr);
    }
    catch (const std::exception& e)
    {
        if (log_normal_enabled())
        {
            std::cerr << "[ADVECTION] numerics vertical advection failed ('"
                      << advection_scheme->name() << "'): " << e.what()
                      << ". Falling back to legacy vertical kernel." << std::endl;
        }
        advect_scalar_1d_z_kernel(src, dst, dt, 0.0);
        return;
    }

    if (tendencies.dqdt_adv.size_r() != NR ||
        tendencies.dqdt_adv.size_th() != NTH ||
        tendencies.dqdt_adv.size_z() != NZ)
    {
        if (log_normal_enabled())
        {
            std::cerr << "[ADVECTION] numerics vertical advection produced unexpected tendency shape. "
                      << "Falling back to legacy vertical kernel." << std::endl;
        }
        advect_scalar_1d_z_kernel(src, dst, dt, 0.0);
        return;
    }

    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NR - 1; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            for (int k = 1; k < NZ - 1; ++k)
            {
                const float q_old = static_cast<float>(src[i][j][k]);
                float dqdt = static_cast<float>(tendencies.dqdt_adv[i][j][k]);
                if (!std::isfinite(static_cast<double>(dqdt)))
                {
                    dqdt = 0.0f;
                }
                float q_new = q_old + static_cast<float>(dt) * dqdt;
                if (!std::isfinite(static_cast<double>(q_new)))
                {
                    q_new = q_old;
                }
                dst[i][j][k] = q_new;
            }
        }
    }
}

/**
 * @brief Advects a scalar field in 3D using directional splitting.
 */
void advect_scalar_3d(Field3D& scalar, double dt, double kappa)
{
    log_runtime_advection_path_once();

    const bool perf_on = perf_enabled();
    const auto scalar_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};

    static Field3D scratch_a;
    static Field3D scratch_b;
    ensure_field_shape(scratch_a);
    ensure_field_shape(scratch_b);

    const double dt_half = dt * 0.5;

    const auto r1_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    advect_scalar_1d_r_kernel(scalar, scratch_a, dt_half, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.r_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - r1_t0).count();
    }

    const auto th1_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    advect_scalar_1d_theta_kernel(scratch_a, scratch_b, dt_half, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.theta_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - th1_t0).count();
    }

    const auto z_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    if (use_numerics_vertical_advection())
    {
        advect_scalar_1d_z_numerics_kernel(scratch_b, scratch_a, dt, kappa);
    }
    else
    {
        advect_scalar_1d_z_kernel(scratch_b, scratch_a, dt, kappa);
    }
    if (perf_on)
    {
        g_advection_perf_totals.z_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - z_t0).count();
    }

    const auto th2_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    advect_scalar_1d_theta_kernel(scratch_a, scratch_b, dt_half, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.theta_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - th2_t0).count();
    }

    const auto r2_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    advect_scalar_1d_r_kernel(scratch_b, scratch_a, dt_half, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.r_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - r2_t0).count();
    }

    const auto diff_t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    apply_diffusion_kernel(scratch_a, scalar, dt, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.diffusion_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - diff_t0).count();
        g_advection_perf_totals.scalar_total_s +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - scalar_t0).count();
        ++g_advection_perf_totals.scalar_calls;
    }
}

/**
 * @brief Advects thermodynamics fields.
 */
void advect_thermodynamics_3d(double dt_advect, double kappa_theta, double kappa_moisture)
{
    const bool perf_on = perf_enabled();
    const auto t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    static int active_nr = -1;
    static int active_nth = -1;
    static int active_nz = -1;
    static bool qc_active = false;
    static bool qr_active = false;
    static bool qi_active = false;
    static bool qs_active = false;
    static bool qh_active = false;
    static bool qg_active = false;

    if (active_nr != NR || active_nth != NTH || active_nz != NZ)
    {
        active_nr = NR;
        active_nth = NTH;
        active_nz = NZ;
        qc_active = false;
        qr_active = false;
        qi_active = false;
        qs_active = false;
        qh_active = false;
        qg_active = false;
    }

    advect_scalar_3d(theta, dt_advect, kappa_theta);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            for (int k = 0; k < NZ; ++k)
            {
                theta(i, j, k) = clamp_theta_k(theta(i, j, k));
            }
        }
    }

    advect_scalar_3d(qv, dt_advect, kappa_moisture);
    qc_active = qc_active || field_has_signal(qc);
    qr_active = qr_active || field_has_signal(qr);
    qi_active = qi_active || field_has_signal(qi);
    qs_active = qs_active || field_has_signal(qs);
    qh_active = qh_active || field_has_signal(qh);
    qg_active = qg_active || field_has_signal(qg);

    if (qc_active) advect_scalar_3d(qc, dt_advect, kappa_moisture);
    if (qr_active) advect_scalar_3d(qr, dt_advect, kappa_moisture);
    if (qi_active) advect_scalar_3d(qi, dt_advect, kappa_moisture);
    if (qs_active) advect_scalar_3d(qs, dt_advect, kappa_moisture);
    if (qh_active) advect_scalar_3d(qh, dt_advect, kappa_moisture);
    if (qg_active) advect_scalar_3d(qg, dt_advect, kappa_moisture);

    int qv_sanitized = 0;
    int hydrometeor_sanitized = 0;
    float* const qv_data = qv.data();
    float* const qc_data = qc_active ? qc.data() : nullptr;
    float* const qr_data = qr_active ? qr.data() : nullptr;
    float* const qi_data = qi_active ? qi.data() : nullptr;
    float* const qs_data = qs_active ? qs.data() : nullptr;
    float* const qh_data = qh_active ? qh.data() : nullptr;
    float* const qg_data = qg_active ? qg.data() : nullptr;
    const std::size_t point_count = qv.size();

    #pragma omp parallel for reduction(+:qv_sanitized, hydrometeor_sanitized)
    for (long long idx = 0; idx < static_cast<long long>(point_count); ++idx)
    {
        const float qv_old = qv_data[idx];
        const float qv_new = clamp_qv_kgkg(qv_old);
        if (qv_new != qv_old)
        {
            ++qv_sanitized;
        }
        qv_data[idx] = qv_new;

        if (qc_data)
        {
            const float old_val = qc_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qc_data[idx] = new_val;
        }
        if (qr_data)
        {
            const float old_val = qr_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qr_data[idx] = new_val;
        }
        if (qi_data)
        {
            const float old_val = qi_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qi_data[idx] = new_val;
        }
        if (qs_data)
        {
            const float old_val = qs_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qs_data[idx] = new_val;
        }
        if (qh_data)
        {
            const float old_val = qh_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qh_data[idx] = new_val;
        }
        if (qg_data)
        {
            const float old_val = qg_data[idx];
            const float new_val = clamp_hydrometeor_kgkg(old_val);
            if (new_val != old_val) { ++hydrometeor_sanitized; }
            qg_data[idx] = new_val;
        }
    }

    if (log_debug_enabled() && (qv_sanitized > 0 || hydrometeor_sanitized > 0))
    {
        std::cout << "[ADVECTION GUARD] qv_sanitized=" << qv_sanitized
                  << ", hydrometeor_sanitized=" << hydrometeor_sanitized << std::endl;
    }

    if (perf_on)
    {
        g_advection_perf_totals.thermo_total_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        ++g_advection_perf_totals.thermo_calls;
    }
}

/**
 * @brief Advects the tracer field.
 */
void advect_tracer_3d(double dt_advect, double kappa)
{
    const bool perf_on = perf_enabled();
    const auto t0 = perf_on ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    advect_scalar_3d(tracer, dt_advect, kappa);
    if (perf_on)
    {
        g_advection_perf_totals.tracer_total_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        ++g_advection_perf_totals.tracer_calls;
    }
}

/**
 * @brief Clears accumulated advection performance statistics.
 */
void reset_advection_perf_stats()
{
    g_advection_perf_totals = AdvectionPerfTotals{};
}

/**
 * @brief Logs aggregated advection performance statistics.
 */
void log_advection_perf_summary()
{
    if (!global_perf_timing_enabled)
    {
        return;
    }

    const double scalar_total = std::max(g_advection_perf_totals.scalar_total_s, 1e-9);
    std::cout << "[ADVECTION PERF] scalar_calls=" << g_advection_perf_totals.scalar_calls
              << ", thermo_calls=" << g_advection_perf_totals.thermo_calls
              << ", tracer_calls=" << g_advection_perf_totals.tracer_calls
              << ", scalar_total_s=" << g_advection_perf_totals.scalar_total_s << std::endl;
    std::cout << "  r_s=" << g_advection_perf_totals.r_s
              << " (" << (100.0 * g_advection_perf_totals.r_s / scalar_total) << "%)" << std::endl;
    std::cout << "  theta_s=" << g_advection_perf_totals.theta_s
              << " (" << (100.0 * g_advection_perf_totals.theta_s / scalar_total) << "%)" << std::endl;
    std::cout << "  z_s=" << g_advection_perf_totals.z_s
              << " (" << (100.0 * g_advection_perf_totals.z_s / scalar_total) << "%)" << std::endl;
    std::cout << "  diffusion_s=" << g_advection_perf_totals.diffusion_s
              << " (" << (100.0 * g_advection_perf_totals.diffusion_s / scalar_total) << "%)" << std::endl;
    std::cout << "  thermo_total_s=" << g_advection_perf_totals.thermo_total_s << std::endl;
    std::cout << "  tracer_total_s=" << g_advection_perf_totals.tracer_total_s << std::endl;
}
