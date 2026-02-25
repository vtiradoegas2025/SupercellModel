/**
 * @file smagorinsky.cpp
 * @brief Implementation for the turbulence module.
 *
 * Provides executable logic for the turbulence runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/turbulence subsystem.
 */

#include "smagorinsky.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>


/**
 * @brief Initializes the Smagorinsky turbulence scheme.
 */
SmagorinskyScheme::SmagorinskyScheme()
    : Cs_(turbulence_constants::C_s_default),
      Pr_t_(turbulence_constants::Pr_t_default),
      Sc_t_(turbulence_constants::Sc_t_default),
      nu_t_max_(1000.0),
      K_max_(1000.0) {
}

/**
 * @brief Gets the required fields.
 */
int SmagorinskyScheme::required_fields() const 
{
    return static_cast<int>(TurbulenceRequirements::BASIC);
}

/**
 * @brief Initializes the Smagorinsky turbulence scheme.
 */
void SmagorinskyScheme::initialize(const TurbulenceConfig& cfg) 
{
    Cs_ = cfg.Cs;
    Pr_t_ = cfg.Pr_t;
    Sc_t_ = cfg.Sc_t;
    nu_t_max_ = cfg.nu_t_max;
    K_max_ = cfg.K_max;

    std::cout << "Initialized Smagorinsky Turbulence:" << std::endl;
    std::cout << "  Cs = " << Cs_ << std::endl;
    std::cout << "  Pr_t = " << Pr_t_ << ", Sc_t = " << Sc_t_ << std::endl;
    std::cout << "  Max nu_t = " << nu_t_max_ << " m²/s" << std::endl;
}

/**
 * @brief Computes the Smagorinsky turbulence scheme.
 */
void SmagorinskyScheme::compute(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state,
    TurbulenceTendencies& tend,
    TurbulenceDiagnostics* diag_opt
) 
{
    tend.dudt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dwdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dthetadt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dqvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dtkedt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);

    compute_eddy_coefficients(cfg, grid, state);

    eddy_viscosity::compute_momentum_diffusion_tendencies(
        cfg, state, grid, nu_t_, tend.dudt_sgs, tend.dvdt_sgs, tend.dwdt_sgs
    );

    eddy_viscosity::compute_scalar_diffusion_tendencies(
        cfg, state, grid, K_theta_, *state.theta, tend.dthetadt_sgs
    );

    
    if (state.qv) 
    {
        eddy_viscosity::compute_scalar_diffusion_tendencies(
            cfg, state, grid, K_q_, *state.qv, tend.dqvdt_sgs
        );
    }

    if (diag_opt) 
    {
        diag_opt->nu_t = nu_t_;
        diag_opt->K_theta = K_theta_;
        diag_opt->K_q = K_q_;
        diag_opt->K_tke.resize(state.NR, state.NTH, state.NZ, 0.0f);
        diag_opt->Cs_eff.resize(state.NR, state.NTH, state.NZ, static_cast<float>(Cs_));
    }
}

/**
 * @brief Computes the eddy coefficients.
 */
void SmagorinskyScheme::compute_eddy_coefficients(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state
) 
{
    nu_t_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_theta_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_q_.resize(state.NR, state.NTH, state.NZ, 0.0f);

    for (int i = 0; i < state.NR; ++i)
    {
        for (int j = 0; j < state.NTH; ++j) 
        {
            for (int k = 0; k < state.NZ; ++k) 
            {
                eddy_viscosity::StrainRate S = eddy_viscosity::compute_strain_rate_3d(state, grid, i, j, k);
                double strain_mag = std::isfinite(S.magnitude) ? S.magnitude : 0.0;

                double Delta = eddy_viscosity::compute_filter_width(cfg, grid, i, j, k);
                if (!std::isfinite(Delta) || Delta <= 0.0)
                {
                    Delta = 1.0;
                }

                double stability_factor = apply_stability_correction(cfg, grid, state, i, j, k);
                if (!std::isfinite(stability_factor) || stability_factor < 0.0)
                {
                    stability_factor = 1.0;
                }

                double nu_t = eddy_viscosity::compute_smagorinsky_viscosity(
                    Cs_, Delta, strain_mag, stability_factor
                );
                if (!std::isfinite(nu_t) || nu_t < 0.0)
                {
                    nu_t = 0.0;
                }

                nu_t = std::min(nu_t, nu_t_max_);

                nu_t_[i][j][k] = static_cast<float>(nu_t);

                double K_theta, K_q, K_tke;
                eddy_viscosity::compute_eddy_diffusivities(nu_t, Pr_t_, Sc_t_, K_theta, K_q, K_tke);
                if (!std::isfinite(K_theta) || K_theta < 0.0) K_theta = 0.0;
                if (!std::isfinite(K_q) || K_q < 0.0) K_q = 0.0;

                K_theta = std::min(K_theta, K_max_);
                K_q = std::min(K_q, K_max_);

                K_theta_[i][j][k] = static_cast<float>(K_theta);
                K_q_[i][j][k] = static_cast<float>(K_q);
            }
        }
    }
}

/**
 * @brief Applies the stability correction.
 */
double SmagorinskyScheme::apply_stability_correction(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state,
    int i, int j, int k
) 
{
    if (cfg.stability_correction == "ri") 
    {
        double N = eddy_viscosity::compute_brunt_vaisala_frequency(state, grid, i, j, k);
        double S_mag = eddy_viscosity::compute_strain_rate_3d(state, grid, i, j, k).magnitude;
        if (!std::isfinite(N) || N < 0.0) N = 0.0;
        if (!std::isfinite(S_mag) || S_mag < 0.0) S_mag = 0.0;

        double Ri = (S_mag > 1e-10) ? (N * N) / (S_mag * S_mag) : 0.0;
        if (!std::isfinite(Ri))
        {
            Ri = 0.0;
        }

        return eddy_viscosity::stability_correction_ri(Ri);
    }
    else 
    {
        return 1.0;
    }
}
