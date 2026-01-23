#include "tornado.hpp"
#include "simulation.hpp"
#include <cmath>
#include <algorithm>

/*This file contains the implementation of the tornado scheme.
This file contains the implementation of the initialize, compute_momentum_tendencies, 
compute_angular_momentum, compute_vorticity_diagnostics, compute_dr, 
compute_dz, and compute_radial_mass_flux functions.*/

TornadoScheme::TornadoScheme()
    : NR_(NR), NTH_(NTH), NZ_(NZ), dr_(dr), dtheta_(dtheta), dz_(dz)
{
}

/*This function computes the momentum tendencies for the tornado scheme.
Takes in the u_r, u_theta, u_z, rho, p, theta, dt, du_r_dt, du_theta_dt, du_z_dt, drho_dt, and dp_dt
and computes the momentum tendencies for the tornado scheme.*/
void TornadoScheme::compute_momentum_tendencies(
    const std::vector<std::vector<std::vector<float>>>& u_r,
    const std::vector<std::vector<std::vector<float>>>& u_theta,
    const std::vector<std::vector<std::vector<float>>>& u_z,
    const std::vector<std::vector<std::vector<float>>>& rho,
    const std::vector<std::vector<std::vector<float>>>& p,
    const std::vector<std::vector<std::vector<float>>>& theta,
    double dt,
    std::vector<std::vector<std::vector<float>>>& du_r_dt,
    std::vector<std::vector<std::vector<float>>>& du_theta_dt,
    std::vector<std::vector<std::vector<float>>>& du_z_dt,
    std::vector<std::vector<std::vector<float>>>& drho_dt,
    std::vector<std::vector<std::vector<float>>>& dp_dt)
{
    // Initialize tendencies to zero at all grid points
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ_; ++k) 
            {
                du_r_dt[i][j][k] = 0.0f;
                du_theta_dt[i][j][k] = 0.0f;
                du_z_dt[i][j][k] = 0.0f;
                drho_dt[i][j][k] = 0.0f;
                dp_dt[i][j][k] = 0.0f;
            }
        }
    }

    // Axisymmetric dynamics - assume independence of θ for mean flow
    // Use j=0 for all calculations (axisymmetric assumption)
    int j = 0;

    // Iterate over all grid points and compute momentum tendencies
    for (int i = 1; i < NR_ - 1; ++i) 
    {
        // Avoid division by zero by adding a small epsilon
        double r = i * dr_ + dynamics_constants::eps;

        // Iterate over all vertical levels
        for (int k = 1; k < NZ_ - 1; ++k) 
        {
            // Get current values (axisymmetric - use j=0)
            double ur = u_r[i][j][k];
            double uth = u_theta[i][j][k];
            double uz = u_z[i][j][k];
            double rho_val = rho[i][j][k];
            double p_val = p[i][j][k];

            // Compute derivatives (axisymmetric)
            double dur_dr = compute_dr(u_r, i, j, k);
            double dur_dz = compute_dz(u_r, i, j, k);

            double duth_dr = compute_dr(u_theta, i, j, k);
            double duth_dz = compute_dz(u_theta, i, j, k);

            double duz_dr = compute_dr(u_z, i, j, k);
            double duz_dz = compute_dz(u_z, i, j, k);

            double dp_dr = compute_dr(p, i, j, k);
            double dp_dz = compute_dz(p, i, j, k);

            // === Radial Momentum Equation ===
            // ∂u_r/∂t = -u_r ∂u_r/∂r - w ∂u_r/∂z + (u_θ²/r) - (1/ρ)∂p/∂r
            double advective_r = -ur * dur_dr - uz * dur_dz;
            double centrifugal = uth * uth / r;
            double pressure_grad_r = -dp_dr / rho_val;

            du_r_dt[i][j][k] = advective_r + centrifugal + pressure_grad_r;

            // Iterate and copy the radial momentum tendency to all azimuthal angles.
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                du_r_dt[i][jj][k] = du_r_dt[i][j][k];
            }

            // === Azimuthal Momentum Equation (tangential wind evolution) ===
            // ∂V/∂t = -U ∂V/∂r - W ∂V/∂z - (U V)/r + F_V
            // where F_V parameterizes turbulent diffusion (negative near V_max)
            double advective_th = -ur * duth_dr - uz * duth_dz;
            double coriolis_th = -ur * uth / r;

            // Simplified turbulent diffusion (negative near max V)
            double fv = 0.0;

            // Check if this is near the radius of maximum wind
            if (i > 0 && i < NR_-1) 
            {
                // Check if this is near the radius of maximum wind
                double v_here = uth;
                double v_inner = (i > 0) ? u_theta[i-1][j][k] : 0.0;
                double v_outer = (i < NR_-1) ? u_theta[i+1][j][k] : 0.0;

                // If the wind speed is greater than the inner and outer wind speeds, apply diffusion.
                if (v_here > v_inner && v_here > v_outer) 
                {
                    // Near V_max, apply diffusion
                    fv = -0.01 * v_here; // negative tendency to reduce max winds
                }
            }

            du_theta_dt[i][j][k] = advective_th + coriolis_th + fv;

            // Iterate and copy the azimuthal momentum tendency to all azimuthal angles.
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                du_theta_dt[i][jj][k] = du_theta_dt[i][j][k];
            }

            // === Vertical Momentum Equation ===
            // ∂w/∂t = -U ∂w/∂r - W ∂w/∂z - (1/ρ)∂p/∂z - g + buoyancy
            double advective_z = -ur * duz_dr - uz * duz_dz;
            double pressure_grad_z = -dp_dz / rho_val;

            // Simplified buoyancy
            double theta_prime = theta[i][j][k] - theta0;
            double buoyancy = dynamics_constants::g * (theta_prime / theta0);

            du_z_dt[i][j][k] = advective_z + pressure_grad_z - dynamics_constants::g + buoyancy;

            // Iterate and copy the vertical momentum tendency to all azimuthal angles.
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                du_z_dt[i][jj][k] = du_z_dt[i][j][k];
            }

            // === Continuity Equation (axisymmetric) ===
            // (1/r) ∂(r u_r)/∂r + ∂w/∂z = 0 (incompressible assumption)
            double drho_dt_val = -rho_val * (dur_dr + ur / r + duz_dz);
            drho_dt[i][j][k] = drho_dt_val;

            // Iterate and copy the density tendency to all azimuthal angles.
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                drho_dt[i][jj][k] = drho_dt_val;
            }

            // === Pressure (simplified for axisymmetric) ===
            // Use cyclostrophic balance near core: V²/r ≈ (1/ρ) ∂p/∂r
            double cyclostrophic_dp_dr = rho_val * uth * uth / r;
            dp_dt[i][j][k] = -ur * dp_dr - uz * dp_dz + cyclostrophic_dp_dr;

            // Iterate and copy the pressure tendency to all azimuthal angles.
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                dp_dt[i][jj][k] = dp_dt[i][j][k];
            }
        }
    }
}

/*This function computes the angular momentum for the tornado scheme.
Takes in the u_r, u_theta, angular_momentum, and angular_momentum_tendency
and computes the angular momentum for the tornado scheme.*/
void TornadoScheme::compute_angular_momentum(
    const std::vector<std::vector<std::vector<float>>>& u_r,
    const std::vector<std::vector<std::vector<float>>>& u_theta,
    std::vector<std::vector<std::vector<float>>>& angular_momentum,
    std::vector<std::vector<std::vector<float>>>& angular_momentum_tendency)
{
    // Iterate over all grid points and compute angular momentum
    for (int i = 0; i < NR_; ++i) 
    {
        // Avoid division by zero by adding a small epsilon
        double r = i * dr_ + dynamics_constants::eps;

        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ_; ++k) 
            {
                double v_theta = u_theta[i][j][k];
                angular_momentum[i][j][k] = r * v_theta;

                // Angular momentum tendency (simplified, ignoring diffusion)
                // DM/Dt = -U ∂M/∂r - W ∂M/∂z
                double ur = u_r[i][j][k];
                double uz = w[i][j][k];

                double dm_dr = 0.0, dm_dz = 0.0;

                // If the row index is greater than 0 and less than the number of rows minus 1, compute the derivative in the radial direction.
                if (i > 0 && i < NR_-1) dm_dr = compute_dr(angular_momentum, i, j, k);
                if (k > 0 && k < NZ_-1) dm_dz = compute_dz(angular_momentum, i, j, k);

                angular_momentum_tendency[i][j][k] = -ur * dm_dr - uz * dm_dz;
            }
        }
    }
}

/*This function computes the vorticity diagnostics for the tornado scheme.
Takes in the u_r, u_theta, u_z, rho, p, vorticity_r, vorticity_theta, vorticity_z, 
stretching_term, tilting_term, and baroclinic_term
and computes the vorticity diagnostics for the tornado scheme.*/
void TornadoScheme::compute_vorticity_diagnostics(
    const std::vector<std::vector<std::vector<float>>>& u_r,
    const std::vector<std::vector<std::vector<float>>>& u_theta,
    const std::vector<std::vector<std::vector<float>>>& u_z,
    const std::vector<std::vector<std::vector<float>>>& rho,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& vorticity_r,
    std::vector<std::vector<std::vector<float>>>& vorticity_theta,
    std::vector<std::vector<std::vector<float>>>& vorticity_z,
    std::vector<std::vector<std::vector<float>>>& stretching_term,
    std::vector<std::vector<std::vector<float>>>& tilting_term,
    std::vector<std::vector<std::vector<float>>>& baroclinic_term)
{
    int j = 0; // axisymmetric

    // Iterate over all grid points and compute vorticity diagnostics
    for (int i = 1; i < NR_ - 1; ++i) 
    {
        // Avoid division by zero by adding a small epsilon
        double r = i * dr_ + dynamics_constants::eps;

        // Iterate over all vertical levels
        for (int k = 1; k < NZ_ - 1; ++k) 
        {
            // Axisymmetric vorticity components
            double dur_dz = compute_dz(u_r, i, j, k);
            double duz_dr = compute_dr(u_z, i, j, k);

            vorticity_r[i][j][k] = 0.0; // ∂u_z/∂θ - ∂u_θ/∂z = 0 (axisymmetric)
            vorticity_theta[i][j][k] = dur_dz - duz_dr; // ∂u_r/∂z - ∂u_z/∂r
            vorticity_z[i][j][k] = (1.0 / r) * compute_dr(u_theta, i, j, k) +
                                   compute_dr(u_theta, i, j, k) / r; // (1/r) ∂(r V)/∂r

            // Iterate over all azimuthal angles and copy to all azimuthal angles
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                vorticity_r[i][jj][k] = vorticity_r[i][j][k];
                vorticity_theta[i][jj][k] = vorticity_theta[i][j][k];
                vorticity_z[i][jj][k] = vorticity_z[i][j][k];
            }

            // Vertical vorticity budget terms (dominant for tornado intensification)
            double zeta = vorticity_z[i][j][k];
            double dw_dz = compute_dz(u_z, i, j, k);

            stretching_term[i][j][k] = zeta * dw_dz;

            // Simplified tilting (axisymmetric approximation)
            tilting_term[i][j][k] = 0.0; // ∂w/∂r * ∂v/∂z - ∂w/∂θ * ∂u/∂z = 0 in axisymmetric

            // Baroclinic generation (simplified)
            double drho_dr = compute_dr(rho, i, j, k);
            double dp_dr = compute_dr(p, i, j, k);
            double rho_sq = rho[i][j][k] * rho[i][j][k];

            // If the density squared is greater than a small epsilon, compute the baroclinic term.
            if (rho_sq > dynamics_constants::eps) 
            {
                baroclinic_term[i][j][k] = (1.0 / rho_sq) * drho_dr * dp_dr;
            } 
            else 
            {
                baroclinic_term[i][j][k] = 0.0;
            }

            // Iterate over all azimuthal angles and copy to all azimuthal angles
            for (int jj = 1; jj < NTH_; ++jj) 
            {
                stretching_term[i][jj][k] = stretching_term[i][j][k];
                tilting_term[i][jj][k] = tilting_term[i][j][k];
                baroclinic_term[i][jj][k] = baroclinic_term[i][j][k];
            }
        }
    }
}

// Helper function implementations



/*This function computes the derivative in the radial direction.
Takes in the field, the row index, the column index, and the level index
and computes the derivative in the radial direction.*/
double TornadoScheme::compute_dr(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const 
{
    return (field[i + 1][j][k] - field[i - 1][j][k]) / (2 * dr_);
}

/*This function computes the derivative in the vertical direction.
Takes in the field, the row index, the column index, and the level index
and computes the derivative in the vertical direction.*/
double TornadoScheme::compute_dz(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const 
{
    return (field[i][j][k + 1] - field[i][j][k - 1]) / (2 * dz_);
}

/*This function computes the radial mass flux.
Takes in the u_r, rho, the row index, and the level index
and computes the radial mass flux.*/
double TornadoScheme::compute_radial_mass_flux(const std::vector<std::vector<std::vector<float>>>& u_r,
                                               const std::vector<std::vector<std::vector<float>>>& rho,
                                               int i, int k) const 
{
    int j = 0; // axisymmetric
    double r = i * dr_ + dynamics_constants::eps;
    return rho[i][j][k] * r * u_r[i][j][k];
}
