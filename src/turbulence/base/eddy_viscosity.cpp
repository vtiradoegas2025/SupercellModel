#include "eddy_viscosity.hpp"
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

/*This file contains the implementation of the eddy viscosity.
It manages the computation of the eddy viscosity.*/

namespace eddy_viscosity 
{

/*This function computes the strain rate.
Takes in the state, grid, and the row, column, and level and computes the strain rate.*/
StrainRate compute_strain_rate_3d(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    int i, int j, int k
) 
{
    StrainRate S;

    // Get grid spacings
    double dx = grid.dx;
    double dy = grid.dy;
    double dz_k = (k >= 0 && k < static_cast<int>(grid.dz.size())) ? grid.dz[k] : grid.dz[0];

    // Velocity gradients (centered differences where possible)
    // du/dx
    double dudx;
    if (i > 0 && i < state.NR - 1) {
        dudx = ((*state.u)[i+1][j][k] - (*state.u)[i-1][j][k]) / (2.0 * dx);
    } else if (i == 0) {
        dudx = ((*state.u)[i+1][j][k] - (*state.u)[i][j][k]) / dx;
    } else { // i == state.NR - 1
        dudx = ((*state.u)[i][j][k] - (*state.u)[i-1][j][k]) / dx;
    }

    // du/dy
    double dudy;
    if (j > 0 && j < state.NTH - 1) {
        dudy = ((*state.u)[i][j+1][k] - (*state.u)[i][j-1][k]) / (2.0 * dy);
    } else if (j == 0) {
        dudy = ((*state.u)[i][j+1][k] - (*state.u)[i][j][k]) / dy;
    } else { // j == state.NTH - 1
        dudy = ((*state.u)[i][j][k] - (*state.u)[i][j-1][k]) / dy;
    }

    // dv/dx
    double dvdx;
    if (i > 0 && i < state.NR - 1) {
        dvdx = ((*state.v)[i+1][j][k] - (*state.v)[i-1][j][k]) / (2.0 * dx);
    } else if (i == 0) {
        dvdx = ((*state.v)[i+1][j][k] - (*state.v)[i][j][k]) / dx;
    } else { // i == state.NR - 1
        dvdx = ((*state.v)[i][j][k] - (*state.v)[i-1][j][k]) / dx;
    }

    // dw/dx
    double dwdx;
    if (i > 0 && i < state.NR - 1) {
        dwdx = ((*state.w)[i+1][j][k] - (*state.w)[i-1][j][k]) / (2.0 * dx);
    } else if (i == 0) {
        dwdx = ((*state.w)[i+1][j][k] - (*state.w)[i][j][k]) / dx;
    } else { // i == state.NR - 1
        dwdx = ((*state.w)[i][j][k] - (*state.w)[i-1][j][k]) / dx;
    }

    // Note: dudv was duplicate of dudy, removing it

    // dv/dy
    double dvdy;
    if (j > 0 && j < state.NTH - 1) {
        dvdy = ((*state.v)[i][j+1][k] - (*state.v)[i][j-1][k]) / (2.0 * dy);
    } else if (j == 0) {
        dvdy = ((*state.v)[i][j+1][k] - (*state.v)[i][j][k]) / dy;
    } else { // j == state.NTH - 1
        dvdy = ((*state.v)[i][j][k] - (*state.v)[i][j-1][k]) / dy;
    }

    // dw/dy
    double dwdy;
    if (j > 0 && j < state.NTH - 1) {
        dwdy = ((*state.w)[i][j+1][k] - (*state.w)[i][j-1][k]) / (2.0 * dy);
    } else if (j == 0) {
        dwdy = ((*state.w)[i][j+1][k] - (*state.w)[i][j][k]) / dy;
    } else { // j == state.NTH - 1
        dwdy = ((*state.w)[i][j][k] - (*state.w)[i][j-1][k]) / dy;
    }

    // du/dz
    double dudz;
    if (k > 0 && k < state.NZ - 1) {
        dudz = ((*state.u)[i][j][k+1] - (*state.u)[i][j][k-1]) / (grid.z_int[k+1] - grid.z_int[k-1]);
    } else if (k == 0) {
        dudz = ((*state.u)[i][j][k+1] - (*state.u)[i][j][k]) / dz_k;
    } else { // k == state.NZ - 1
        dudz = ((*state.u)[i][j][k] - (*state.u)[i][j][k-1]) / dz_k;
    }

    // dv/dz
    double dvdz;
    if (k > 0 && k < state.NZ - 1) {
        dvdz = ((*state.v)[i][j][k+1] - (*state.v)[i][j][k-1]) / (grid.z_int[k+1] - grid.z_int[k-1]);
    } else if (k == 0) {
        dvdz = ((*state.v)[i][j][k+1] - (*state.v)[i][j][k]) / dz_k;
    } else { // k == state.NZ - 1
        dvdz = ((*state.v)[i][j][k] - (*state.v)[i][j][k-1]) / dz_k;
    }

    // dw/dz
    double dwdz;
    if (k > 0 && k < state.NZ - 1) {
        dwdz = ((*state.w)[i][j][k+1] - (*state.w)[i][j][k-1]) / (grid.z_int[k+1] - grid.z_int[k-1]);
    } else if (k == 0) {
        dwdz = ((*state.w)[i][j][k+1] - (*state.w)[i][j][k]) / dz_k;
    } else { // k == state.NZ - 1
        dwdz = ((*state.w)[i][j][k] - (*state.w)[i][j][k-1]) / dz_k;
    }

    // Strain rate tensor components (symmetric)
    S.S11 = dudx;           // du/dx
    S.S12 = 0.5 * (dudy + dvdx);  // 0.5*(du/dy + dv/dx)
    S.S13 = 0.5 * (dudz + dwdx);  // 0.5*(du/dz + dw/dx)

    S.S21 = S.S12;          // symmetric
    S.S22 = dvdy;           // dv/dy
    S.S23 = 0.5 * (dvdz + dwdy);  // 0.5*(dv/dz + dw/dy)

    S.S31 = S.S13;          // symmetric
    S.S32 = S.S23;          // symmetric
    S.S33 = dwdz;           // dw/dz

    // Strain rate magnitude: |S| = sqrt(2 * Sij * Sij)
    S.magnitude = std::sqrt(2.0 * (
        S.S11*S.S11 + S.S22*S.S22 + S.S33*S.S33 +
        2.0*(S.S12*S.S12 + S.S13*S.S13 + S.S23*S.S23)
    ));

    return S;
}

/*This function computes the filter width.
Takes in the configuration, grid, and the row, column, and level and computes the filter width.*/
double compute_filter_width(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    int i, int j, int k
) 
{

    // If the filter width is "dx", compute the filter width.
    if (cfg.filter_width == "dx") 
    {
        // Simple horizontal grid spacing
        return std::sqrt(grid.dx * grid.dy);
    }
    
    // If the filter width is "cubic_root", compute the filter width.
    else if (cfg.filter_width == "cubic_root") 
    {
        // LES-style cubic root (default)
        double dz_k = (k >= 0 && k < static_cast<int>(grid.dz.size())) ? grid.dz[k] : grid.dz[0];
        return std::cbrt(grid.dx * grid.dy * dz_k);
    }
    else 
    {
        // Default fallback
        return std::cbrt(grid.dx * grid.dy * 100.0);  // assume 100m vertical spacing
    }
}

/*This function computes the Smagorinsky viscosity.
Takes in the Smagorinsky constant, filter width, strain magnitude, and stability factor and computes the Smagorinsky viscosity.*/

double compute_smagorinsky_viscosity(
    double Cs,
    double Delta,
    double strain_mag,
    double stability_factor
) 
{
    // Smagorinsky formula: nu_t = (Cs * Delta)^2 * |S|
    double nu_t = Cs * Cs * Delta * Delta * strain_mag;

    // Apply stability correction
    nu_t *= stability_factor;

    // Ensure non-negative
    return std::max(nu_t, 0.0);
}

/*This function computes the eddy diffusivities.
Takes in the eddy viscosity, turbulent Prandtl numbers, and the eddy diffusivities and computes the eddy diffusivities.*/
void compute_eddy_diffusivities(
    double nu_t,
    double Pr_t,
    double Sc_t,
    double& K_theta,
    double& K_q,
    double& K_tke
) 
{
    // Eddy diffusivities from eddy viscosity and turbulent Prandtl numbers
    K_theta = nu_t / Pr_t;      // temperature diffusivity
    K_q = nu_t / Sc_t;          // moisture diffusivity
    K_tke = nu_t / Pr_t;        // TKE diffusivity (same as temperature)
}

/*This function computes the scalar diffusion tendency.
Takes in the state, grid, the eddy diffusivities, and the scalar and computes the scalar diffusion tendency.*/
double compute_scalar_diffusion_tendency(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& K_field,
    const Field3D& phi,
    int i, int j, int k,
    int var_index
) 
{
    // This is a simplified implementation - in practice, you'd need proper
    // flux-divergence computation across the entire domain

    // For now, return zero (implement proper diffusion in full scheme)
    return 0.0;
}

/*This function computes the momentum diffusion tendencies.
Takes in the state, grid, the eddy viscosity, and the momentum and computes the momentum diffusion tendencies.*/
void compute_momentum_diffusion_tendencies(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& nu_t,
    Field3D& dudt_sgs,
    Field3D& dvdt_sgs,
    Field3D& dwdt_sgs
) 
{
    // Simplified momentum diffusion (would need proper implementation)
    // For now, apply weak damping to prevent instability

    // Iterate over the rows, columns, and levels and compute the momentum diffusion tendencies.        
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < state.NR; ++i) 
    {
        // Iterate over the columns and compute the momentum diffusion tendencies.
        for (int j = 0; j < state.NTH; ++j)
        {
            // Iterate over the levels and compute the momentum diffusion tendencies.
            for (int k = 0; k < state.NZ; ++k) 
            {
                // Compute the local eddy viscosity.
                double nu_local = static_cast<double>(nu_t[i][j][k]);
                double rho_local = static_cast<double>((*state.rho)[i][j][k]);

                // Simple damping: -nu * u / rho (simplified diffusion)
                dudt_sgs[i][j][k] = -nu_local * static_cast<double>((*state.u)[i][j][k]) / rho_local;
                dvdt_sgs[i][j][k] = -nu_local * static_cast<double>((*state.v)[i][j][k]) / rho_local;
                dwdt_sgs[i][j][k] = -nu_local * static_cast<double>((*state.w)[i][j][k]) / rho_local;
            }
        }
    }
}

/*This function computes the scalar diffusion tendencies.
Takes in the state, grid, the eddy diffusivities, and the scalar and computes the scalar diffusion tendencies.*/
void compute_scalar_diffusion_tendencies(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& K_field,
    const Field3D& phi,
    Field3D& dphi_dt_sgs
) 
{
    // Simplified scalar diffusion (would need proper implementation)
    // Iterate over the rows, columns, and levels and compute the scalar diffusion tendencies.
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < state.NR; ++i) 
    {
        // Iterate over the columns and compute the scalar diffusion tendencies.
        for (int j = 0; j < state.NTH; ++j) 
        {
            // Iterate over the levels and compute the scalar diffusion tendencies.
            for (int k = 0; k < state.NZ; ++k) 
            {
                // Compute the local eddy diffusivity.
                double K_local = static_cast<double>(K_field[i][j][k]);
                double phi_local = static_cast<double>(phi[i][j][k]);

                // Simple damping: -K * phi (simplified diffusion)
                dphi_dt_sgs[i][j][k] = -K_local * phi_local;
            }
        }
    }
}

double stability_correction_ri(double Ri, double Ri_crit) {
    // Reduce mixing in stable stratification
    if (Ri > 0.0) {
        return std::max(0.1, 1.0 - Ri / Ri_crit);
    } else {
        // Slightly enhance mixing in unstable conditions
        return 1.0 - 0.5 * Ri;
    }
}

/*This function computes the TKE mixing length.
Takes in the filter width, TKE, Brunt-Väisälä frequency, and the stability factor and computes the TKE mixing length.*/
double compute_tke_mixing_length(
    double Delta,
    double e,
    double N,
    double c_l
) 
{
    // Stability-limited mixing length
    double l_grid = Delta;  // grid-based limit
    double l_stability = (N > 1e-6) ? c_l * std::sqrt(e) / N : l_grid;

    return std::min(l_grid, l_stability);
}

/*This function computes the Brunt-Väisälä frequency.
Takes in the state, and the row, column, and level and computes the Brunt-Väisälä frequency.*/
double compute_brunt_vaisala_frequency(
    const TurbulenceStateView& state,
    int i, int j, int k
)
{
    // If the theta is not set or the level is greater than the number of levels, return zero.
    if (!state.theta || k >= state.NZ - 1) return 0.0;

    // Get the theta at the current level.
    double theta_k = static_cast<double>((*state.theta)[i][j][k]);
    // Get the theta at the next level.
    double theta_kp1 = static_cast<double>((*state.theta)[i][j][k+1]);

    double dtheta_dz = (theta_kp1 - theta_k) / 100.0;  // assume 100m layer

    // Brunt-Väisälä frequency: N² = (g/θ) * dθ/dz
    double N2 = (turbulence_constants::g / theta_k) * dtheta_dz;

    return (N2 > 0.0) ? std::sqrt(N2) : 0.0;
}

/*This function applies the positivity limits.
Takes in the field, minimum value, and maximum value and applies the positivity limits.*/
void apply_positivity_limits(
    Field3D& field,
    double min_value,
    double max_value
) 
{
    // Iterate over all grid points and apply the positivity limits
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < field.size_r(); ++i) 
    {
        for (int j = 0; j < field.size_th(); ++j) 
        {
            for (int k = 0; k < field.size_z(); ++k) 
            {
                field[i][j][k] = std::max(static_cast<float>(min_value),
                      std::min(static_cast<float>(max_value), static_cast<float>(field[i][j][k])));
            }
        }
    }
}

/*This function initializes the 3D field.
Takes in the field, number of rows, number of columns, number of levels, and the value and initializes the 3D field.*/
void initialize_3d_field(
    std::vector<std::vector<std::vector<float>>>& field,
    int NR, int NTH, int NZ,
    float value
) 
{
    // Assign the field.
    field.assign(NR, std::vector<std::vector<float>>(
                NTH, std::vector<float>(NZ, value)));
}

} // namespace eddy_viscosity
