#include "weno5.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This class implements the WENO5 advection scheme.
It is a subclass of the AdvectionSchemeBase class.
It implements the compute_flux_divergence method.*/
WENO5Scheme::WENO5Scheme() : weno_epsilon_(numerics_constants::weno_epsilon) 
{
}

/*This function initializes the WENO5 advection scheme.
It is a subclass of the AdvectionSchemeBase class.
It implements the compute_flux_divergence method.*/
void WENO5Scheme::initialize() 
{
    initialize(AdvectionConfig{});
}

/*This function initializes the WENO5 advection scheme.
It is a subclass of the AdvectionSchemeBase class.
It implements the compute_flux_divergence method.*/
void WENO5Scheme::initialize(const AdvectionConfig& cfg) 
{
    config_ = cfg;
    weno_epsilon_ = numerics_constants::weno_epsilon;

    std::cout << "Initialized WENO5 advection scheme" << std::endl;

    // If the positivity preservation is enabled, print a message.
    if (cfg.positivity) 
    {
        std::cout << "  Positivity preservation enabled" << std::endl;
    }
}

/*This function computes the flux divergence.
Takes in the configuration, state, tendencies, and diagnostics and computes the flux divergence.*/
void WENO5Scheme::compute_flux_divergence(
    const AdvectionConfig& cfg,
    const AdvectionStateView& state,
    AdvectionTendencies& tendencies,
    AdvectionDiagnostics* diag_opt
) 
{
    // If the state is invalid, throw an error.
    if (!state.q || !state.grid) 
    {
        throw std::runtime_error("Invalid state for WENO5 advection");
    }

    const int NR = state.q->size();
    const int NTH = (*state.q)[0].size();
    const int NZ = (*state.q)[0][0].size();

    // Initialize tendencies
    tendencies.dqdt_adv.assign(NR, std::vector<std::vector<double>>(
        NTH, std::vector<double>(NZ, 0.0)));

    double max_cfl = 0.0;

    // Simplified: advect in vertical direction only
    // Full implementation would need directional splitting

    // Iterate over the horizontal columns and compute the flux divergence.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and compute the flux divergence.
        for (int j = 0; j < NTH; ++j) 
        {
            // Extract 1D column data
            std::vector<double> q_col(NZ);
            std::vector<double> w_col(NZ);
            for (int k = 0; k < NZ; ++k) {
                q_col[k] = (*state.q)[i][j][k];
                w_col[k] = (*state.w)[i][j][k];
            }

            // Compute vertical advection with WENO5
            std::vector<double> dqdt_col(NZ, 0.0);
            weno5_advect_1d(q_col, w_col, state.grid->dz[0], 1.0, dqdt_col);

            // Iterate over the vertical levels and store the tendencies.
            for (int k = 0; k < NZ; ++k)
            {
                tendencies.dqdt_adv[i][j][k] = dqdt_col[k];
            }

            // Iterate over the vertical levels and compute the CFL.
            for (int k = 0; k < NZ; ++k) 
            {
                double cfl = std::abs(w_col[k]) / state.grid->dz[k];
                max_cfl = std::max(max_cfl, cfl);
            }
        }
    }

    // If the diagnostics are requested, store the maximum CFL and the suggested time step.
    if (diag_opt) 
    {
        diag_opt->max_cfl_z = max_cfl;
        diag_opt->suggested_dt = cfg.cfl_target * state.grid->dz[0] / max_cfl;
    }
}

/*This function suggests the time step.
Takes in the configuration and state and suggests the time step.*/
double WENO5Scheme::suggest_dt(
    const AdvectionConfig& cfg,
    const AdvectionStateView& state
) 
{
    // If the grid is invalid, return 1.0.
    if (!state.grid) return 1.0;

    const int NR = state.q->size();
    const int NTH = (*state.q)[0].size();
    const int NZ = (*state.q)[0][0].size();

    double max_cfl = 0.0;

    // Iterate over the horizontal columns and compute the CFL in all directions.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and compute the CFL in the x direction.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the vertical levels and compute the CFL in the z direction.
            for (int k = 0; k < NZ; ++k) 
            {
                // If the u velocity is valid, compute the CFL in the x direction.
                if (state.u) 
                {
                    double cfl_x = std::abs((*state.u)[i][j][k]) / state.grid->dx;
                    max_cfl = std::max(max_cfl, cfl_x);
                }

                // If the v velocity is valid, compute the CFL in the y direction.
                if (state.v) 
                {
                    double cfl_y = std::abs((*state.v)[i][j][k]) / state.grid->dy;
                    max_cfl = std::max(max_cfl, cfl_y);
                }

                // If the w velocity is valid, compute the CFL in the z direction.
                if (state.w) 
                {
                    double cfl_z = std::abs((*state.w)[i][j][k]) / state.grid->dz[k];
                    max_cfl = std::max(max_cfl, cfl_z);
                }
            }
        }
    }

    return cfg.cfl_target / max_cfl;
}


/*This function performs the WENO5 reconstruction for the left state.
Takes in the state and the index and performs the WENO5 reconstruction for the left state.*/
double WENO5Scheme::weno5_reconstruct_left(const std::vector<double>& q, int i) 
{
    // If the index is less than 2 or greater than the number of grid points minus 3, return the value at the index.
    if (i < 2 || i >= static_cast<int>(q.size()) - 3) 
    {
        // Fall back to simple reconstruction near boundaries
        return q[i];
    }

    // Candidate reconstructions
    double q0 = q_tilde_0(q, i);
    double q1 = q_tilde_1(q, i);
    double q2 = q_tilde_2(q, i);

    // Smoothness indicators
    double beta0 = smoothness_indicator(q, i-2, i+1);
    double beta1 = smoothness_indicator(q, i-1, i+2);
    double beta2 = smoothness_indicator(q, i, i+3);

    // Nonlinear weights
    double alpha0 = d0 / ((weno_epsilon_ + beta0) * (weno_epsilon_ + beta0));
    double alpha1 = d1 / ((weno_epsilon_ + beta1) * (weno_epsilon_ + beta1));
    double alpha2 = d2 / ((weno_epsilon_ + beta2) * (weno_epsilon_ + beta2));

    double sum_alpha = alpha0 + alpha1 + alpha2;
    double omega0 = alpha0 / sum_alpha;
    double omega1 = alpha1 / sum_alpha;
    double omega2 = alpha2 / sum_alpha;

    // Weighted reconstruction
    return omega0 * q0 + omega1 * q1 + omega2 * q2;
}

/*This function performs the WENO5 reconstruction for the right state.
Takes in the state and the index and performs the WENO5 reconstruction for the right state.*/
double WENO5Scheme::weno5_reconstruct_right(const std::vector<double>& q, int i) 
{
    // Symmetric reconstruction for right state
    // For right-biased reconstruction at i+1/2, use points i-1 to i+4
    if (i < 1 || i >= static_cast<int>(q.size()) - 4) {
        return q[i+1];
    }

    // Shift index for right reconstruction
    double q0 = q_tilde_0(q, i+1);
    double q1 = q_tilde_1(q, i+1);
    double q2 = q_tilde_2(q, i+1);

    double beta0 = smoothness_indicator(q, i-1, i+2);
    double beta1 = smoothness_indicator(q, i, i+3);
    double beta2 = smoothness_indicator(q, i+1, i+4);

    double alpha0 = d0 / ((weno_epsilon_ + beta0) * (weno_epsilon_ + beta0));
    double alpha1 = d1 / ((weno_epsilon_ + beta1) * (weno_epsilon_ + beta1));
    double alpha2 = d2 / ((weno_epsilon_ + beta2) * (weno_epsilon_ + beta2));

    double sum_alpha = alpha0 + alpha1 + alpha2;
    double omega0 = alpha0 / sum_alpha;
    double omega1 = alpha1 / sum_alpha;
    double omega2 = alpha2 / sum_alpha;

    return omega0 * q0 + omega1 * q1 + omega2 * q2;
}

/*This function computes the smoothness indicator.
Takes in the state and the start and end indices and computes the smoothness indicator.*/
double WENO5Scheme::smoothness_indicator(const std::vector<double>& q, int start, int end) 
{
    // If the start is less than 0 or the end is greater than the number of grid points or the end minus the start is not 3, return a large value.
    if (start < 0 || end >= static_cast<int>(q.size()) || end - start != 3) 
    {
        return 1e6; // Large value for invalid stencils
    }

    double beta = 0.0;
    beta += 13.0/12.0 * std::pow(q[start] - 2*q[start+1] + q[start+2], 2);
    beta += 1.0/4.0 * std::pow(q[start] - 4*q[start+1] + 3*q[start+2], 2);

    return beta;
}

/*This function computes the q_tilde_0.
Takes in the state and the index and computes the q_tilde_0.*/
double WENO5Scheme::q_tilde_0(const std::vector<double>& q, int i) 
{
    // Stencil 0: q[i-2], q[i-1], q[i+1]
    return (1.0/3.0)*q[i-2] - (7.0/6.0)*q[i-1] + (11.0/6.0)*q[i];
}

/*This function computes the q_tilde_1.
Takes in the state and the index and computes the q_tilde_1.*/
double WENO5Scheme::q_tilde_1(const std::vector<double>& q, int i) 
{
    // Stencil 1: q[i-1], q[i], q[i+1]
    return -(1.0/6.0)*q[i-1] + (5.0/6.0)*q[i] + (1.0/3.0)*q[i+1];
}

/*This function computes the q_tilde_2.
Takes in the state and the index and computes the q_tilde_2.*/
double WENO5Scheme::q_tilde_2(const std::vector<double>& q, int i) 
{
    // Stencil 2: q[i], q[i+1], q[i+2]
    return (1.0/3.0)*q[i] + (5.0/6.0)*q[i+1] - (1.0/6.0)*q[i+2];
}

/*This function computes the Rusanov flux.
Takes in the left and right states and the velocity and computes the Rusanov flux.*/
double WENO5Scheme::rusanov_flux(double q_left, double q_right, double velocity) 
{
    // Rusanov (local Lax-Friedrichs) flux for advection
    double alpha = std::abs(velocity);
    return 0.5 * velocity * (q_left + q_right) - 0.5 * alpha * (q_right - q_left);
}

/*This function advects the 1D field.
Takes in the state, the velocity, the grid spacing, the time step, 
and the tendencies and advects the 1D field.*/
void WENO5Scheme::weno5_advect_1d(
    const std::vector<double>& q,
    const std::vector<double>& velocity,
    double dx, double dt,
    std::vector<double>& dqdt
) 
{
    const int n = q.size();
    dqdt.assign(n, 0.0);

    // Compute fluxes at interfaces
    std::vector<double> flux(n+1, 0.0);

    // Iterate over the grid points and compute the fluxes.
    for (int i = 0; i < n-1; ++i) 
    {
        // Reconstruct left and right states at interface i+1/2
        double q_left = weno5_reconstruct_left(q, i);
        double q_right = weno5_reconstruct_right(q, i);

        // Compute numerical flux
        flux[i+1] = rusanov_flux(q_left, q_right, velocity[i]);
    }

    // Apply boundary conditions (zero flux)
    flux[0] = 0.0;
    flux[n] = 0.0;

    // Iterate over the grid points and compute the flux divergence.
    for (int i = 0; i < n; ++i) 
    {
        double dflux_dx = (flux[i+1] - flux[i]) / dx;
        dqdt[i] = -dflux_dx;
    }

    // If the positivity preservation is enabled, apply the positivity preservation.
    if (config_.positivity) 
    {
        apply_positivity_preservation(dqdt);
    }
}

/*This function applies the positivity preservation.
Takes in the state and the minimum value and applies the positivity preservation.*/
void WENO5Scheme::apply_positivity_preservation(std::vector<double>& q, double min_value) 
{
    // Iterate over the grid points and apply the positivity preservation.
    for (auto& val : q) 
    {
        // If the value is less than the minimum value, set the value to the minimum value.
        if (val < min_value) 
        {
            val = min_value;
        }
    }
}
