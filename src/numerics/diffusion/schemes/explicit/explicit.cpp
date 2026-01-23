#include "explicit.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This class implements the explicit diffusion scheme.
It is a subclass of the DiffusionSchemeBase class.
It implements the compute_diffusion_tendencies method.*/

/*This function initializes the explicit diffusion scheme.
Takes in the configuration and initializes the explicit diffusion scheme.*/

ExplicitDiffusionScheme::ExplicitDiffusionScheme() 
{
}

/*This function initializes the explicit diffusion scheme.
Takes in the configuration and initializes the explicit diffusion scheme.*/
void ExplicitDiffusionScheme::initialize() 
{
    initialize(DiffusionConfig{});
}

/*This function initializes the explicit diffusion scheme.
Takes in the configuration and initializes the explicit diffusion scheme.*/
void ExplicitDiffusionScheme::initialize(const DiffusionConfig& cfg) 
{
    config_ = cfg;

    std::cout << "Initialized explicit diffusion scheme" << std::endl;
    std::cout << "  K_h = " << cfg.K_h << " m²/s, K_v = " << cfg.K_v << " m²/s" << std::endl;
}

/*This function computes the diffusion tendencies.
Takes in the configuration, state, tendencies, and diagnostics and computes the diffusion tendencies.*/
void ExplicitDiffusionScheme::compute_diffusion_tendencies(
    const DiffusionConfig& cfg,
    const DiffusionStateView& state,
    DiffusionTendencies& tendencies,
    DiffusionDiagnostics* diag_opt
) 
{
    // If the grid is invalid, throw an error.
    if (!state.grid) {throw std::runtime_error("Invalid grid for explicit diffusion");}

    const int NR = state.theta->size();
    const int NTH = (*state.theta)[0].size();
    const int NZ = (*state.theta)[0][0].size();

    // Initialize tendencies
    tendencies.dudt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dvdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dwdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dthetadt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dqvdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));

    // Create diffusivity fields
    std::vector<std::vector<std::vector<double>>> K_scalar(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, cfg.K_h)));
    std::vector<std::vector<std::vector<double>>> K_momentum(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, cfg.K_h)));

    // If the variable diffusivity is enabled and the variable diffusivity is valid, copy the variable diffusivity.
    if (cfg.use_variable_K && state.K_scalar) 
    {
        // Iterate over the horizontal columns and vertical levels and copy the variable diffusivity.
        for (int i = 0; i < NR; ++i) 
        {
            // Iterate over the horizontal columns and vertical levels and copy the variable diffusivity.
            for (int j = 0; j < NTH; ++j) 
            {
                // Iterate over the vertical levels and copy the variable diffusivity.
                for (int k = 0; k < NZ; ++k) 
                {
                    K_scalar[i][j][k] = (*state.K_scalar)[i][j][k];
                    K_momentum[i][j][k] = state.K_momentum ? (*state.K_momentum)[i][j][k] : cfg.K_h;
                }
            }
        }
    }

    // If the u, v, and w velocities are valid, compute the momentum diffusion.
    if (state.u && state.v && state.w) 
    {
        compute_momentum_diffusion(*state.u, *state.v, *state.w, K_momentum,
                                 *state.grid, tendencies.dudt_diff, tendencies.dvdt_diff, tendencies.dwdt_diff);
    }

    // If the theta field is valid, compute the scalar diffusion.
    if (state.theta) 
    {
        compute_scalar_diffusion(*state.theta, K_scalar, *state.grid, tendencies.dthetadt_diff);
    }

    // If the qv field is valid, compute the scalar diffusion.
    if (state.qv) 
    {
        compute_scalar_diffusion(*state.qv, K_scalar, *state.grid, tendencies.dqvdt_diff);
    }

    // If the diagnostics are requested, store the effective diffusivity and the maximum diffusion number.
    if (diag_opt) 
    {
        diag_opt->K_effective = K_scalar;
        diag_opt->max_diffusion_number = check_stability(cfg, state);
    }
}

/*This function checks the stability of the diffusion scheme.
Takes in the configuration and state and checks the stability of the diffusion scheme.*/
double ExplicitDiffusionScheme::check_stability(
    const DiffusionConfig& cfg,
    const DiffusionStateView& state
) 
{
    // If the grid is invalid, return 0.0.
    if (!state.grid) return 0.0;

    double max_diff_num = 0.0;
    double K_max = std::max(cfg.K_h, cfg.K_v);

    // Check diffusion number in all directions
    double diff_num_x = 2.0 * K_max * cfg.dt_diffusion / (state.grid->dx * state.grid->dx);
    double diff_num_y = 2.0 * K_max * cfg.dt_diffusion / (state.grid->dy * state.grid->dy);
    double diff_num_z = 2.0 * K_max * cfg.dt_diffusion / (state.grid->dz[0] * state.grid->dz[0]);

    max_diff_num = std::max({diff_num_x, diff_num_y, diff_num_z});

    return max_diff_num;
}

/*This function computes the scalar diffusion.
Takes in the field, the diffusivity field, the grid, and the tendency and computes the scalar diffusion.*/
void ExplicitDiffusionScheme::compute_scalar_diffusion(
    const std::vector<std::vector<std::vector<double>>>& field,
    const std::vector<std::vector<std::vector<double>>>& K_field,
    const GridMetrics& grid,
    std::vector<std::vector<std::vector<double>>>& tendency
) 
{
    // If the grid is invalid, throw an error.
    const int NR = field.size();
    const int NTH = field[0].size();
    const int NZ = field[0][0].size();

    // Simplified 3D diffusion (vertical only for stability)
    // Full 3D would need proper flux-divergence computation

    // Iterate over the horizontal columns and vertical levels and compute the scalar diffusion.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and vertical levels and compute the scalar diffusion.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the vertical levels and compute the scalar diffusion.
            for (int k = 1; k < NZ - 1; ++k) 
            {
                // Compute the vertical grid spacing.
                double dz = grid.dz[k];
                double dz_up = grid.dz[k-1];
                double dz_down = grid.dz[k+1];

                // Vertical diffusion: d/dz(K dphi/dz)
                double dphi_dz_up = (field[i][j][k] - field[i][j][k-1]) / dz_up;
                double dphi_dz_down = (field[i][j][k+1] - field[i][j][k]) / dz_down;

                double K_avg_up = 0.5 * (K_field[i][j][k] + K_field[i][j][k-1]);
                double K_avg_down = 0.5 * (K_field[i][j][k] + K_field[i][j][k+1]);

                double flux_up = -K_avg_up * dphi_dz_up;
                double flux_down = -K_avg_down * dphi_dz_down;

                double dflux_dz = (flux_down - flux_up) / dz;

                tendency[i][j][k] = dflux_dz;
            }
        }
    }

    // Iterate over the horizontal columns and vertical levels and set the boundary conditions.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and vertical levels and set the boundary conditions.
        for (int j = 0; j < NTH; ++j) 
        {
            // Set the boundary conditions.
            tendency[i][j][0] = 0.0;
            tendency[i][j][NZ-1] = 0.0;
        }
    }
}

/*This function computes the momentum diffusion.
Takes in the u, v, w velocities, the diffusivity field, the grid, 
and the tendency and computes the momentum diffusion.*/
void ExplicitDiffusionScheme::compute_momentum_diffusion(
    const std::vector<std::vector<std::vector<double>>>& u,
    const std::vector<std::vector<std::vector<double>>>& v,
    const std::vector<std::vector<std::vector<double>>>& w,
    const std::vector<std::vector<std::vector<double>>>& nu_t,
    const GridMetrics& grid,
    std::vector<std::vector<std::vector<double>>>& du_dt,
    std::vector<std::vector<std::vector<double>>>& dv_dt,
    std::vector<std::vector<std::vector<double>>>& dw_dt
) 
{
    const int NR = u.size();
    const int NTH = u[0].size();
    const int NZ = u[0][0].size();

    
    // Simplified vertical momentum diffusion

    // Iterate over the horizontal columns and vertical levels and compute the momentum diffusion.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and vertical levels and compute the momentum diffusion.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the vertical levels and compute the momentum diffusion.
            for (int k = 1; k < NZ - 1; ++k) {
                double dz = grid.dz[k];

                // Vertical diffusion for u
                double du_dz_up = (u[i][j][k] - u[i][j][k-1]) / grid.dz[k-1];
                double du_dz_down = (u[i][j][k+1] - u[i][j][k]) / grid.dz[k+1];
                double nu_avg = 0.5 * (nu_t[i][j][k] + nu_t[i][j][k-1]);
                du_dt[i][j][k] = nu_avg * (du_dz_down - du_dz_up) / dz;

                // Vertical diffusion for v
                double dv_dz_up = (v[i][j][k] - v[i][j][k-1]) / grid.dz[k-1];
                double dv_dz_down = (v[i][j][k+1] - v[i][j][k]) / grid.dz[k+1];
                dv_dt[i][j][k] = nu_avg * (dv_dz_down - dv_dz_up) / dz;

                // Vertical diffusion for w
                double dw_dz_up = (w[i][j][k] - w[i][j][k-1]) / grid.dz[k-1];
                double dw_dz_down = (w[i][j][k+1] - w[i][j][k]) / grid.dz[k+1];
                dw_dt[i][j][k] = nu_avg * (dw_dz_down - dw_dz_up) / dz;
            }
        }
    }

    // Iterate over the horizontal columns and vertical levels and set the boundary conditions.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and vertical levels and set the boundary conditions.
        for (int j = 0; j < NTH; ++j) 
        {
            du_dt[i][j][0] = 0.0;
            dv_dt[i][j][0] = 0.0;
            dw_dt[i][j][0] = 0.0;
            du_dt[i][j][NZ-1] = 0.0;
            dv_dt[i][j][NZ-1] = 0.0;
            dw_dt[i][j][NZ-1] = 0.0;
        }
    }
}
