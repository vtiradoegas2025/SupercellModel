#include "tvd.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This constructor initializes the TVD scheme with a default limiter function.*/
TVDScheme::TVDScheme() : limiter_function_(nullptr) 
{
}

/*This function initializes the TVD scheme with a default configuration.*/
void TVDScheme::initialize() 
{
    initialize(AdvectionConfig{});
}

/*This function initializes the TVD scheme with a configuration.
Takes in the configuration and initializes the TVD scheme with the configuration.*/
void TVDScheme::initialize(const AdvectionConfig& cfg) 
{
    // Set the configuration.
    config_ = cfg;

    // If the limiter id is minmod, set the limiter function to the minmod limiter function.
    if (cfg.limiter_id == "minmod") 
    {
        limiter_function_ = &minmod_limiter;
    } 

    // If the limiter id is vanleer, set the limiter function to the vanleer limiter function.
    else if (cfg.limiter_id == "vanleer") 
    {
        limiter_function_ = &vanleer_limiter;
    } 
    
    // If the limiter id is superbee, set the limiter function to the superbee limiter function.
    else if (cfg.limiter_id == "superbee") 
    {
        limiter_function_ = &superbee_limiter;
    } 

    // If the limiter id is mc, set the limiter function to the mc limiter function.
    else if (cfg.limiter_id == "mc") 
    {
        limiter_function_ = &mc_limiter;
    } 
    
    // If the limiter id is universal, set the limiter function to the universal limiter function.
    else if (cfg.limiter_id == "universal") {
        limiter_function_ = &universal_limiter;
    } 
    else 
    {
        std::cerr << "Warning: Unknown limiter '" << cfg.limiter_id << "', using MC limiter" << std::endl;
        limiter_function_ = &mc_limiter;
    }

    std::cout << "Initialized TVD advection scheme with " << cfg.limiter_id << " limiter" << std::endl;
}

// TVD limiter functions

/*This function computes the minmod limiter.
Takes in the slope ratio and computes the minmod limiter.*/
double TVDScheme::minmod_limiter(double r) 
{
    return std::max(0.0, std::min(1.0, r));
}


/*This function computes the vanleer limiter.
Takes in the slope ratio and computes the vanleer limiter.*/
double TVDScheme::vanleer_limiter(double r) 
{
    return (std::abs(r) + r) / (1.0 + std::abs(r));
}


/*This function computes the superbee limiter.
Takes in the slope ratio and computes the superbee limiter.*/
double TVDScheme::superbee_limiter(double r) 
{
    return std::max({0.0, std::min(1.0, 2.0 * r), std::min(2.0, r)});
}


/*This function computes the mc limiter.
Takes in the slope ratio and computes the mc limiter.*/
double TVDScheme::mc_limiter(double r) 
{
    return std::max(0.0, std::min({(1.0 + r) / 2.0, 2.0, 2.0 * r}));
}


/*This function computes the universal limiter.
Takes in the slope ratio and computes the universal limiter.*/
double TVDScheme::universal_limiter(double r) 
{
    // Universal limiter (shape-preserving)
    double phi = 0.0;

    // If the slope ratio is between 0 and 1, compute the universal limiter.
    if (r >= 0.0 && r <= 1.0) 
    {
        phi = std::min(2.0 * r, (1.0 + r) / 2.0);
    } 

    // If the slope ratio is greater than 1, compute the universal limiter.
    else if (r > 1.0) 
    {
        phi = std::min(r, 2.0);
    }
    return phi;
}


/*This function computes the flux divergence.
Takes in the configuration, state, tendencies, and diagnostics and computes the flux divergence.*/
void TVDScheme::compute_flux_divergence(
    const AdvectionConfig& cfg,
    const AdvectionStateView& state,
    AdvectionTendencies& tendencies,
    AdvectionDiagnostics* diag_opt
) 
{
    // If the state is invalid, throw an error.
    if (!state.q || !state.grid) {throw std::runtime_error("Invalid state for TVD advection");}

    const int NR = state.q->size();
    const int NTH = (*state.q)[0].size();
    const int NZ = (*state.q)[0][0].size();

    // Initialize tendencies
    tendencies.dqdt_adv.assign(NR, std::vector<std::vector<double>>(
        NTH, std::vector<double>(NZ, 0.0)));

    double max_cfl = 0.0;

    // For each horizontal column, advect in vertical direction
    // (This is a simplified implementation - full 3D would need directional splitting)
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NTH; ++j) {
            // Extract 1D column data
            std::vector<double> q_col(NZ);
            std::vector<double> w_col(NZ);
            for (int k = 0; k < NZ; ++k) {
                q_col[k] = (*state.q)[i][j][k];
                w_col[k] = (*state.w)[i][j][k];
            }

            // Compute vertical advection
            std::vector<double> dqdt_col(NZ, 0.0);
            advect_1d(q_col, w_col, state.grid->dz[0], 1.0, dqdt_col);  // dt=1.0 for tendency

            // Iterate over the vertical levels and store the tendencies.
            for (int k = 0; k < NZ; ++k) 
            {
                tendencies.dqdt_adv[i][j][k] = dqdt_col[k];
            }

            //Iterate over the vertical levels and compute the CFL.
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
double TVDScheme::suggest_dt(
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


/*This function performs the MUSCL reconstruction.
Takes in the state, and the limiter function and performs the MUSCL reconstruction.*/
void TVDScheme::muscl_reconstruct(
    const std::vector<double>& q,
    std::vector<double>& q_left,
    std::vector<double>& q_right,
    double (*limiter)(double)
) 
{
    // Get the number of grid points.
    const int n = q.size();

    // Resize the left and right states.
    q_left.resize(n);
    q_right.resize(n);

    // Iterate over the grid points and perform the MUSCL reconstruction.
    for (int i = 1; i < n - 1; ++i) 
    {
        // Compute slope ratio r
        double denominator = q[i+1] - q[i] + numerics_constants::epsilon;
        double r = (q[i] - q[i-1]) / denominator;

        // Apply limiter
        double phi = limiter(r);

        // Compute limited slope
        double delta_q = phi * (q[i+1] - q[i]);

        // Reconstruct left and right states
        q_left[i] = q[i] - 0.5 * delta_q;
        q_right[i] = q[i] + 0.5 * delta_q;
    }

    // Boundary conditions (first-order upwind)
    q_left[0] = q[0];
    q_right[0] = q[0];
    q_left[n-1] = q[n-1];
    q_right[n-1] = q[n-1];
}


/*This function computes the numerical flux.
Takes in the left and right states and the velocity and computes the numerical flux.*/
double TVDScheme::numerical_flux(double q_left, double q_right, double velocity) 
{
    // If the velocity is positive, compute the upwind flux.
    if (velocity >= 0.0) 
    {
        return velocity * q_left;
    } 
    else 
    {
        return velocity * q_right;
    }
}

/*This function advects the 1D field.
Takes in the state, the velocity, the grid spacing, the time step, 
and the tendencies and advects the 1D field.*/
void TVDScheme::advect_1d(
    const std::vector<double>& q,
    const std::vector<double>& velocity,
    double dx, double dt,
    std::vector<double>& dqdt
) 
{
    const int n = q.size();
    dqdt.assign(n, 0.0);

    // MUSCL reconstruction
    std::vector<double> q_left, q_right;
    muscl_reconstruct(q, q_left, q_right, limiter_function_);

    // Iterate over the grid points and compute the fluxes and tendencies.
    for (int i = 0; i < n - 1; ++i) 
    {
        // Flux at interface i+1/2
        double flux_right = numerical_flux(q_right[i], q_left[i+1], velocity[i]);
        double flux_left = (i > 0) ? numerical_flux(q_right[i-1], q_left[i], velocity[i-1]) : 0.0;

        // Flux divergence
        double dflux_dx = (flux_right - flux_left) / dx;

        // Tendency
        dqdt[i] -= dflux_dx;
    }

    // If the positivity limiter is enabled, apply the positivity limiter.
    if (config_.positivity) 
    {
        apply_positivity_limiter(dqdt);
    }
}


/*This function applies the positivity limiter.
Takes in the state and the minimum value and applies the positivity limiter.*/
void TVDScheme::apply_positivity_limiter(std::vector<double>& q, double min_value) 
{
    // Iterate over the grid points and apply the positivity limiter.
    for (auto& val : q) 
    {
        val = std::max(val, min_value);
    }
}
