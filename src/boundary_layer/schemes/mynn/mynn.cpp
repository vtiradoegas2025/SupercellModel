/**
 * @file mynn.cpp
 * @brief Implementation for the boundary_layer module.
 *
 * Provides executable logic for the boundary_layer runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/boundary_layer subsystem.
 */

#include "mynn.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>


/**
 * @brief Initializes the MYNN boundary layer scheme.
 */
MYNNScheme::MYNNScheme()
    : a1_(1.18), a2_(0.665), b1_(24.0), b2_(15.0),
      c1_(0.137), c2_(0.75), c3_(0.352), c4_(0.0), c5_(0.2),
      ce1_(1.0), ce2_(1.33), lmax_(500.0) {
}

/**
 * @brief Returns the required fields for the MYNN boundary layer scheme.
 */
int MYNNScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE) |
           static_cast<int>(BoundaryLayerRequirements::TKE);
}

/**
 * @brief Initializes the MYNN boundary layer scheme Takes in the boundary layer configuration.
 */
void MYNNScheme::initialize(const BoundaryLayerConfig& cfg) 
{
    lmax_ = cfg.pbl_max_height * 0.5;

    std::cout << "Initialized MYNN Boundary Layer:" << std::endl;
    std::cout << "  Prognostic TKE: enabled" << std::endl;
    std::cout << "  Max mixing length: " << lmax_ << " m" << std::endl;
}

/**
 * @brief Computes the column for the MYNN boundary layer scheme.
 */
void MYNNScheme::compute_column(
    const BoundaryLayerConfig& cfg,
    const SurfaceConfig& sfc,
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    BoundaryLayerDiagnostics* diag_opt
) 
{
    const size_t nz = col.theta->size();
    if (nz == 0 || !col.u || !col.v || !col.theta || !col.qv || !col.p || !col.rho || !col.z_int || col.z_int->size() < nz + 1)
    {
        return;
    }

    std::vector<double> tke_col(nz, 0.1);
    if (col.tke && col.tke->size() == nz)
    {
        for (size_t k = 0; k < nz; ++k)
        {
            const double tke_k = (*col.tke)[k];
            tke_col[k] = std::isfinite(tke_k) ? std::max(tke_k, 0.001) : 0.1;
        }
    }

    tend.dudt_pbl.assign(nz, 0.0);
    tend.dvdt_pbl.assign(nz, 0.0);
    tend.dthetadt_pbl.assign(nz, 0.0);
    tend.dqvdt_pbl.assign(nz, 0.0);
    tend.dtkedt_pbl.assign(nz, 0.0);

    surface_fluxes::BulkFluxes fluxes;

    if (cfg.apply_surface_fluxes) 
    {
        fluxes = surface_fluxes::compute_surface_fluxes(
            cfg,
            sfc, col.u_sfc, col.v_sfc, col.theta_sfc, col.qv_sfc,
            (*col.theta)[0], (*col.qv)[0], (*col.p)[0], (*col.rho)[0], col.z_sfc
        );
    }

    double h = diagnose_pbl_height(col, tke_col, cfg);

    std::vector<double> l_mix(nz);
    compute_mixing_length(col, h, l_mix);

    std::vector<double> K_m(nz + 1, 0.0);
    std::vector<double> K_h(nz + 1, 0.0);
    std::vector<double> K_e(nz + 1, 0.0);
    compute_eddy_diffusivities(col, tke_col, l_mix, K_m, K_h, K_e);

    update_tke(col, tke_col, l_mix, fluxes, K_m, K_h, K_e, cfg.dt_pbl, tend.dtkedt_pbl);

    apply_diffusion_tendencies(col, K_m, K_h, tend);

    if (diag_opt) 
    {
        diag_opt->pbl_height = h;
        diag_opt->ustar = fluxes.ustar;
        diag_opt->tau_u = fluxes.tau_u;
        diag_opt->tau_v = fluxes.tau_v;
        diag_opt->H = fluxes.H;
        diag_opt->E = fluxes.E;
        diag_opt->Cd = fluxes.Cd;
        diag_opt->Ch = fluxes.Ch;
        diag_opt->Ce = fluxes.Ce;
        diag_opt->K_m = K_m;
        diag_opt->K_h = K_h;
        diag_opt->K_e = K_e;
    }
}

/**
 * @brief Diagnoses the PBL height.
 */
double MYNNScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& tke_col,
    const BoundaryLayerConfig& cfg
) 
{
    const double min_h = 100.0;
    const double max_h = cfg.pbl_max_height;

    const double tke_min = 0.01;

    for (size_t k = 1; k < col.theta->size(); ++k)
    {
        if (tke_col[k] < tke_min) 
        {
            double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
            return std::max(min_h, std::min(max_h, z_level));
        }
    }

    return max_h;
}

/**
 * @brief Computes the mixing length.
 */
void MYNNScheme::compute_mixing_length(
    const BoundaryLayerColumnStateView& col,
    double h,
    std::vector<double>& l_mix
) 
{
    const size_t nz = col.theta->size();
    const double kappa = boundary_layer_constants::kappa;
    const double h_safe = std::max(h, 1.0);

    for (size_t k = 0; k < nz; ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k + 1]);

        double l_blackadar = kappa * z_level / (1.0 + kappa * z_level / lmax_);

        double l_wall = kappa * z_level;
        double l_top = kappa * std::max(h_safe - z_level, 0.0);

        l_mix[k] = std::min({l_blackadar, l_wall, l_top, lmax_});
        if (!std::isfinite(l_mix[k]))
        {
            l_mix[k] = 1.0;
        }
        l_mix[k] = std::max(l_mix[k], 1.0);
    }
}

/**
 * @brief Computes the eddy diffusivities.
 */
void MYNNScheme::compute_eddy_diffusivities(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& tke_col,
    const std::vector<double>& l_mix,
    std::vector<double>& K_m,
    std::vector<double>& K_h,
    std::vector<double>& K_e
) {
    const size_t nz = col.theta->size();

    for (size_t k = 1; k <= nz; ++k) 
    {
        size_t k_cell = k - 1;

        double e_k = std::max(tke_col[k_cell], 1.0e-8);

        double q2 = 2.0 * e_k;
        double gh = 0.0;

        const double denom = (1.0 - gh) * (1.0 + a2_ * gh);
        double sm = 0.01;
        double sh = 0.01;
        if (std::isfinite(denom) && std::abs(denom) > 1.0e-6)
        {
            sm = a1_ * (1.0 - c1_ * gh) / denom;
            sh = a1_ * (1.0 - c2_ * gh) / denom;
        }

        sm = std::max(sm, 0.01);
        sh = std::max(sh, 0.01);

        double l_k = std::max(l_mix[k_cell], 1.0);
        const double sqrt_q2 = std::sqrt(std::max(q2, 0.0));
        K_m[k] = l_k * sqrt_q2 * sm;
        K_h[k] = l_k * sqrt_q2 * sh;
        K_e[k] = l_k * sqrt_q2 * sh;

        if (!std::isfinite(K_m[k])) K_m[k] = 0.0;
        if (!std::isfinite(K_h[k])) K_h[k] = 0.0;
        if (!std::isfinite(K_e[k])) K_e[k] = 0.0;
        K_m[k] = std::clamp(K_m[k], 0.0, 1000.0);
        K_h[k] = std::clamp(K_h[k], 0.0, 1000.0);
        K_e[k] = std::clamp(K_e[k], 0.0, 1000.0);
    }

    K_m[0] = 0.0;
    K_h[0] = 0.0;
    K_e[0] = 0.0;
}

/**
 * @brief Updates the TKE field.
 */
void MYNNScheme::update_tke(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& tke_col,
    const std::vector<double>& l_mix,
    const surface_fluxes::BulkFluxes& fluxes,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    const std::vector<double>& K_e,
    double dt,
    std::vector<double>& tke_tend
) 
{
    const size_t nz = col.theta->size();
    const double dt_safe = std::max(dt, 1.0e-6);
    constexpr double max_abs_tke_step = 5.0;
    const double max_abs_tke_tendency = max_abs_tke_step / dt_safe;


    for (size_t k = 0; k < nz; ++k) 
    {
        double dz_k = std::max((*col.z_int)[k+1] - (*col.z_int)[k], 1.0e-6);
        double du_dz = ((*col.u)[k] - (k > 0 ? (*col.u)[k-1] : 0.0)) / dz_k;
        double dv_dz = ((*col.v)[k] - (k > 0 ? (*col.v)[k-1] : 0.0)) / dz_k;
        double S2 = du_dz * du_dz + dv_dz * dv_dz;

        double P_s = K_m[k+1] * S2;

        double dtheta_dz = ((*col.theta)[k] - (k > 0 ? (*col.theta)[k-1] : (*col.theta)[k])) / dz_k;
        const double theta_safe = std::max(std::abs((*col.theta)[k]), 1.0e-6);
        double P_b = -K_h[k+1] * dtheta_dz * boundary_layer_constants::g / theta_safe;

        double q2 = 2.0 * std::max(tke_col[k], 1.0e-8);
        double q3 = q2 * sqrt(q2);
        const double l_mix_k = std::max(l_mix[k], 1.0);
        double eps = q3 / (b1_ * l_mix_k);

        double flux_up = 0.0;
        if (k < nz - 1)
        {
            const double dz_up = std::max((*col.z_int)[k + 2] - (*col.z_int)[k + 1], 1.0e-6);
            const double de_dz_up = (tke_col[k + 1] - tke_col[k]) / dz_up;
            flux_up = -std::max(K_e[k + 1], 0.0) * de_dz_up;
        }
        double flux_dn = 0.0;
        if (k > 0)
        {
            const double dz_dn = std::max((*col.z_int)[k + 1] - (*col.z_int)[k], 1.0e-6);
            const double de_dz_dn = (tke_col[k] - tke_col[k - 1]) / dz_dn;
            flux_dn = -std::max(K_e[k], 0.0) * de_dz_dn;
        }
        const double diffusion = -(flux_up - flux_dn) / dz_k;
        const double rho_safe = std::max((*col.rho)[k], 1.0e-6);

        double tke_tendency = P_s + P_b - eps + diffusion / rho_safe;

        if (k == 0)
        {
            const double z_sfc_safe = std::max(col.z_sfc, 1.0);
            const double ustar = std::max(fluxes.ustar, 0.0);
            double ustar3 = ustar * ustar * ustar;
            tke_tendency += ustar3 / (boundary_layer_constants::kappa * z_sfc_safe);
        }

        if (!std::isfinite(tke_tendency))
        {
            tke_tendency = 0.0;
        }
        tke_tend[k] = std::clamp(tke_tendency, -max_abs_tke_tendency, max_abs_tke_tendency);
    }
}

/**
 * @brief Applies the diffusion tendencies.
 */
void MYNNScheme::apply_diffusion_tendencies(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    BoundaryLayerTendencies& tend
) 
{
    const size_t nz = col.theta->size();
    constexpr double max_abs_momentum_tendency = 20.0;
    constexpr double max_abs_theta_tendency = 5.0;
    constexpr double max_abs_qv_tendency = 0.01;

    for (size_t k = 0; k < nz; ++k) 
    {
        double dz_k = std::max((*col.z_int)[k+1] - (*col.z_int)[k], 1.0e-6);
        const double rho_safe = std::max((*col.rho)[k], 1.0e-6);

        double K_m_avg = 0.5 * (K_m[k] + K_m[k+1]);
        double du_dz = (k < nz-1) ?
            ((*col.u)[k+1] - (*col.u)[k]) / dz_k : 0.0;
        double dv_dz = (k < nz-1) ?
            ((*col.v)[k+1] - (*col.v)[k]) / dz_k : 0.0;

        tend.dudt_pbl[k] = K_m_avg * du_dz / dz_k / rho_safe;
        tend.dvdt_pbl[k] = K_m_avg * dv_dz / dz_k / rho_safe;

        double K_h_avg = 0.5 * (K_h[k] + K_h[k+1]);
        double dtheta_dz = (k < nz-1) ?
            ((*col.theta)[k+1] - (*col.theta)[k]) / dz_k : 0.0;
        double dqv_dz = (k < nz-1) ?
            ((*col.qv)[k+1] - (*col.qv)[k]) / dz_k : 0.0;

        tend.dthetadt_pbl[k] = K_h_avg * dtheta_dz / dz_k / rho_safe;
        tend.dqvdt_pbl[k] = K_h_avg * dqv_dz / dz_k / rho_safe;

        if (!std::isfinite(tend.dudt_pbl[k])) tend.dudt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dvdt_pbl[k])) tend.dvdt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dthetadt_pbl[k])) tend.dthetadt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dqvdt_pbl[k])) tend.dqvdt_pbl[k] = 0.0;

        tend.dudt_pbl[k] = std::clamp(tend.dudt_pbl[k], -max_abs_momentum_tendency, max_abs_momentum_tendency);
        tend.dvdt_pbl[k] = std::clamp(tend.dvdt_pbl[k], -max_abs_momentum_tendency, max_abs_momentum_tendency);
        tend.dthetadt_pbl[k] = std::clamp(tend.dthetadt_pbl[k], -max_abs_theta_tendency, max_abs_theta_tendency);
        tend.dqvdt_pbl[k] = std::clamp(tend.dqvdt_pbl[k], -max_abs_qv_tendency, max_abs_qv_tendency);
    }
}
