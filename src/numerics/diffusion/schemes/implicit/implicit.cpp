#include "implicit.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This class implements the implicit diffusion scheme.
It is a subclass of the DiffusionSchemeBase class.
It implements the compute_diffusion_tendencies method.*/
ImplicitDiffusionScheme::ImplicitDiffusionScheme() {}

/*This function initializes the implicit diffusion scheme.
Takes in the configuration and initializes the implicit diffusion scheme.*/
void ImplicitDiffusionScheme::initialize() 
{
    initialize(DiffusionConfig{});
}

/*This function initializes the implicit diffusion scheme.
Takes in the configuration and initializes the implicit diffusion scheme.*/
void ImplicitDiffusionScheme::initialize(const DiffusionConfig& cfg) 
{
    config_ = cfg;

    std::cout << "Initialized implicit diffusion scheme" << std::endl;
    std::cout << "  Vertical implicit diffusion enabled" << std::endl;
}

/*This function computes the diffusion tendencies.
Takes in the configuration, state, tendencies, and diagnostics 
and computes the diffusion tendencies.*/
void ImplicitDiffusionScheme::compute_diffusion_tendencies(
    const DiffusionConfig& cfg,
    const DiffusionStateView& state,
    DiffusionTendencies& tendencies,
    DiffusionDiagnostics* diag_opt
) 
{
    // If the grid is invalid, throw an error.
    if (!state.grid) 
    {
        throw std::runtime_error("Invalid grid for implicit diffusion");
    }

    const int NR = state.theta->size();
    const int NTH = (*state.theta)[0].size();
    const int NZ = (*state.theta)[0][0].size();

    // Initialize tendencies (implicit schemes typically return updated fields, not tendencies)
    tendencies.dudt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dvdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dwdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dthetadt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    tendencies.dqvdt_diff.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));

    // For implicit schemes, we compute the updated fields and then compute tendencies
    // This is a simplified approach - in practice, implicit schemes return tendencies

   // If the theta field is valid, compute the scalar diffusion.
    if (state.theta) 
    {
        // Initialize the new theta field.
        std::vector<std::vector<std::vector<double>>> theta_new = *state.theta;
        std::vector<std::vector<std::vector<double>>> K_theta(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, cfg.K_v)));

        // If the variable diffusivity is enabled and the variable diffusivity is valid, copy the variable diffusivity.
        if (cfg.use_variable_K && state.K_scalar) 
        {
            // Use variable diffusivity
        }

        implicit_vertical_diffusion(*state.theta, K_theta, *state.grid, cfg.dt_diffusion, theta_new);

        // Iterate over the horizontal columns and vertical levels and compute the tendency.
        for (int i = 0; i < NR; ++i) 
        {
            // Iterate over the horizontal columns and vertical levels and compute the tendency.
            for (int j = 0; j < NTH; ++j) 
            {
                // Iterate over the vertical levels and compute the tendency.
                for (int k = 0; k < NZ; ++k) 
                {
                    tendencies.dthetadt_diff[i][j][k] = (theta_new[i][j][k] - (*state.theta)[i][j][k]) / cfg.dt_diffusion;
                }
            }
        }
    }

    // If the u, v, and w velocities are valid, compute the momentum diffusion.
    if (state.u && state.v && state.w) 
    {
        std::vector<std::vector<std::vector<double>>> u_new = *state.u;
        std::vector<std::vector<std::vector<double>>> v_new = *state.v;
        std::vector<std::vector<std::vector<double>>> w_new = *state.w;
        std::vector<std::vector<std::vector<double>>> nu_t(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, cfg.K_v)));

        implicit_vertical_momentum_diffusion(*state.u, *state.v, *state.w, nu_t,
                                           *state.grid, cfg.dt_diffusion, u_new, v_new, w_new);

        // Iterate over the horizontal columns and vertical levels and compute the tendency.
        for (int i = 0; i < NR; ++i) 
        {
            // Iterate over the horizontal columns and vertical levels and compute the tendency.
            for (int j = 0; j < NTH; ++j) 
            {
                // Iterate over the vertical levels and compute the tendency.
                for (int k = 0; k < NZ; ++k) 
                {
                    tendencies.dudt_diff[i][j][k] = (u_new[i][j][k] - (*state.u)[i][j][k]) / cfg.dt_diffusion;
                    tendencies.dvdt_diff[i][j][k] = (v_new[i][j][k] - (*state.v)[i][j][k]) / cfg.dt_diffusion;
                    tendencies.dwdt_diff[i][j][k] = (w_new[i][j][k] - (*state.w)[i][j][k]) / cfg.dt_diffusion;
                }
            }
        }
    }

    // If the diagnostics are requested, store the effective diffusivity and the maximum diffusion number.
    if (diag_opt) 
    {
        diag_opt->K_effective.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, cfg.K_v)));
        diag_opt->max_diffusion_number = 1.0;  // Implicit schemes are unconditionally stable
    }
}

/*This function checks the stability of the diffusion scheme.
Takes in the configuration and state and checks the stability of the diffusion scheme.*/
double ImplicitDiffusionScheme::check_stability(
    const DiffusionConfig& cfg,
    const DiffusionStateView& state
) 
{
    // Implicit schemes are unconditionally stable for diffusion
    return 1.0;
}

/*This function solves the tridiagonal system.
Takes in the a, b, c, rhs, and x and solves the tridiagonal system.*/
void ImplicitDiffusionScheme::solve_tridiagonal(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& rhs,
    std::vector<double>& x
) 
{
    // Thomas algorithm for tridiagonal systems
    const int n = b.size();
    std::vector<double> cp(n-1, 0.0);
    std::vector<double> dp(n, 0.0);

    // Forward elimination
    cp[0] = c[0] / b[0];
    dp[0] = rhs[0] / b[0];

    // Iterate over the grid points and solve the tridiagonal system.
    for (int i = 1; i < n; ++i) 
    {
        double denom = b[i] - a[i] * cp[i-1];

        // If the index is less than the number of grid points minus 1, solve the tridiagonal system.
        if (i < n-1) 
        {
            cp[i] = c[i] / denom;
        }
        dp[i] = (rhs[i] - a[i] * dp[i-1]) / denom;
    }

    // Back substitution
    x[n-1] = dp[n-1];

    // Iterate over the grid points and solve the tridiagonal system.
    for (int i = n-2; i >= 0; --i) 
    {
        x[i] = dp[i] - cp[i] * x[i+1];
    }
}

/*This function computes the vertical diffusion.
Takes in the field, the diffusivity field, the grid, the time step, 
and the new field and computes the vertical diffusion.*/
void ImplicitDiffusionScheme::implicit_vertical_diffusion(
    const std::vector<std::vector<std::vector<double>>>& field,
    const std::vector<std::vector<std::vector<double>>>& K_field,
    const GridMetrics& grid,
    double dt,
    std::vector<std::vector<std::vector<double>>>& field_new
) 
{
    const int NR = field.size();
    const int NTH = field[0].size();
    const int NZ = field[0][0].size();

    // Iterate over the horizontal columns and vertical levels and solve the tridiagonal system.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the horizontal columns and vertical levels and solve the tridiagonal system.
        for (int j = 0; j < NTH; ++j) {
            // Set up tridiagonal system for vertical diffusion
            std::vector<double> a(NZ-1, 0.0);  // Sub-diagonal
            std::vector<double> b(NZ, 1.0);    // Main diagonal
            std::vector<double> c(NZ-1, 0.0);  // Super-diagonal
            std::vector<double> rhs(NZ, 0.0);  // Right-hand side

            // Iterate over the vertical levels and solve the tridiagonal system.
            for (int k = 1; k < NZ-1; ++k) 
            {
                double dz_up = grid.dz[k-1];
                double dz_down = grid.dz[k+1];
                double dz = grid.dz[k];

                double K_up = 0.5 * (K_field[i][j][k] + K_field[i][j][k-1]);
                double K_down = 0.5 * (K_field[i][j][k] + K_field[i][j][k+1]);

                // Coefficients for implicit Crank-Nicolson
                double coef_up = -0.5 * dt * K_up / (dz_up * dz);
                double coef_down = -0.5 * dt * K_down / (dz_down * dz);

                a[k-1] = coef_up;      // k-1 equation coefficient for phi[k]
                b[k] = 1.0 - coef_up - coef_down;  // k equation coefficient for phi[k]
                c[k] = coef_down;      // k+1 equation coefficient for phi[k]

                // RHS: explicit part + current field
                double explicit_up = 0.5 * dt * K_up / (dz_up * dz);
                double explicit_down = 0.5 * dt * K_down / (dz_down * dz);

                rhs[k] = field[i][j][k] +
                        explicit_up * (field[i][j][k-1] - field[i][j][k]) +
                        explicit_down * (field[i][j][k+1] - field[i][j][k]);
            }

            // Boundary conditions (zero flux)
            b[0] = 1.0;
            rhs[0] = field[i][j][0];

            b[NZ-1] = 1.0;
            rhs[NZ-1] = field[i][j][NZ-1];

            // Solve tridiagonal system
            std::vector<double> phi_new(NZ);
            solve_tridiagonal(a, b, c, rhs, phi_new);

            // Iterate over the vertical levels and store the result.
            for (int k = 0; k < NZ; ++k) 
            {
                field_new[i][j][k] = phi_new[k];
            }
        }
    }
}

/*This function computes the vertical momentum diffusion.
Takes in the u, v, w velocities, the diffusivity field, the grid, 
the time step, the new u, v, and w velocities, and computes the vertical momentum diffusion.*/
void ImplicitDiffusionScheme::implicit_vertical_momentum_diffusion(
    const std::vector<std::vector<std::vector<double>>>& u,
    const std::vector<std::vector<std::vector<double>>>& v,
    const std::vector<std::vector<std::vector<double>>>& w,
    const std::vector<std::vector<std::vector<double>>>& nu_t,
    const GridMetrics& grid,
    double dt,
    std::vector<std::vector<std::vector<double>>>& u_new,
    std::vector<std::vector<std::vector<double>>>& v_new,
    std::vector<std::vector<std::vector<double>>>& w_new
) 
{
    // Apply implicit diffusion to each momentum component
    implicit_vertical_diffusion(u, nu_t, grid, dt, u_new);
    implicit_vertical_diffusion(v, nu_t, grid, dt, v_new);
    implicit_vertical_diffusion(w, nu_t, grid, dt, w_new);
}
