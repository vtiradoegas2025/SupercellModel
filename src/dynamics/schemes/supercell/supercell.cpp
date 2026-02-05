#include "supercell.hpp"
#include "simulation.hpp"
#include <cmath>
#include <algorithm>

/*This function initializes the supercell scheme.
Takes in the number of rows, the number of azimuthal angles, the number of vertical levels,
the radial grid spacing, the azimuthal grid spacing, and the vertical grid spacing
and initializes the supercell scheme.*/
SupercellScheme::SupercellScheme()
    : NR_(NR), NTH_(NTH), NZ_(NZ),
      dr_(dr), dtheta_(dtheta), dz_(dz) 
{
    // Constructor initializes grid parameters from global simulation parameters
}

/*This function computes the momentum tendencies for the supercell scheme.
Takes in the u_r, u_theta, u_z, rho, p, theta, dt, du_r_dt, du_theta_dt, du_z_dt, drho_dt, and dp_dt
and computes the momentum tendencies for the supercell scheme.*/

void SupercellScheme::compute_momentum_tendencies(
    const Field3D& u_r,
    const Field3D& u_theta,
    const Field3D& u_z,
    const Field3D& rho,
    const Field3D& p,
    const Field3D& theta,
    double dt,
    Field3D& du_r_dt,
    Field3D& du_theta_dt,
    Field3D& du_z_dt,
    Field3D& drho_dt,
    Field3D& dp_dt)
{
    // Iterate over all grid points and initialize tendencies to zero
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

    // Iterate over all grid points and compute tendencies
    for (int i = 1; i < NR_ - 1; ++i) 
    {
        // Avoid division by zero by adding a small epsilon
        double r = i * dr_ + dynamics_constants::eps; 

        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH_; ++j) 
        {
            // Get the previous and next azimuthal angles.
            int j_prev = (j - 1 + NTH_) % NTH_;
            int j_next = (j + 1) % NTH_;

            // Iterate over all vertical levels
            for (int k = 1; k < NZ_ - 1; ++k) 
            {
                // Get current values of u_r, u_theta, u_z, rho, p, and theta.
                double ur = u_r[i][j][k];
                double uth = u_theta[i][j][k];
                double uz = u_z[i][j][k];
                double rho_val = rho[i][j][k];
                double p_val = p[i][j][k];
                double theta_val = theta[i][j][k];

                // Compute derivatives
                double dur_dr = compute_dr(u_r, i, j, k);
                double dur_dth = compute_dtheta(u_r, i, j, k);
                double dur_dz = compute_dz(u_r, i, j, k);

                double duth_dr = compute_dr(u_theta, i, j, k);
                double duth_dth = compute_dtheta(u_theta, i, j, k);
                double duth_dz = compute_dz(u_theta, i, j, k);

                double duz_dr = compute_dr(u_z, i, j, k);
                double duz_dth = compute_dtheta(u_z, i, j, k);
                double duz_dz = compute_dz(u_z, i, j, k);

                double dp_dr = compute_dr(p, i, j, k);
                double dp_dth = compute_dtheta(p, i, j, k);
                double dp_dz = compute_dz(p, i, j, k);

                // === Radial Momentum Equation ===
                // Du_r/Dt = -u·∇u_r + (u_θ²/r) - (1/ρ)∂p/∂r + F_r
                double advective_r = -ur * dur_dr - (uth / r) * dur_dth - uz * dur_dz;
                double centrifugal = uth * uth / r;
                double pressure_grad_r = -dp_dr / rho_val;

                du_r_dt[i][j][k] = advective_r + centrifugal + pressure_grad_r;

                // === Azimuthal Momentum Equation ===
                // Du_θ/Dt = -u·∇u_θ - (u_r u_θ)/r - (1/(ρ r))∂p/∂θ + F_θ
                double advective_th = -ur * duth_dr - (uth / r) * duth_dth - uz * duth_dz;
                double coriolis_th = -ur * uth / r;
                double pressure_grad_th = -dp_dth / (rho_val * r);

                du_theta_dt[i][j][k] = advective_th + coriolis_th + pressure_grad_th;

                // === Vertical Momentum Equation ===
                // Du_z/Dt = -u·∇u_z - (1/ρ)∂p/∂z - g + buoyancy + F_z
                double advective_z = -ur * duz_dr - (uth / r) * duz_dth - uz * duz_dz;
                double pressure_grad_z = -dp_dz / rho_val;

                // Buoyancy: g * (θ'/θ̄) where θ' ≈ θ - θ0
                double theta_prime = theta_val - theta0;
                double buoyancy = dynamics_constants::g * (theta_prime / theta0);

                du_z_dt[i][j][k] = advective_z + pressure_grad_z - dynamics_constants::g + buoyancy;

                // === Mass Continuity (compressible) ===
                // Dρ/Dt = -ρ ∇·u
                double divergence = dur_dr + duz_dz + ur / r + duth_dth / r;
                drho_dt[i][j][k] = -rho_val * divergence;

                // === Pressure (ideal gas law + energy conservation) ===
                // Simplified: pressure tendency based on thermodynamic energy
                double gamma_term = dynamics_constants::gamma * p_val * divergence;
                double advection_p = -ur * dp_dr - (uth / r) * dp_dth - uz * dp_dz;
                dp_dt[i][j][k] = -gamma_term + advection_p;
            }
        }
    }
}

/*This function computes the vorticity diagnostics for the supercell scheme.
Takes in the u_r, u_theta, u_z, rho, p, vorticity_r, vorticity_theta, vorticity_z, stretching_term, tilting_term, and baroclinic_term
and computes the vorticity diagnostics for the supercell scheme.*/
void SupercellScheme::compute_vorticity_diagnostics(
    const Field3D& u_r,
    const Field3D& u_theta,
    const Field3D& u_z,
    const Field3D& rho,
    const Field3D& p,
    Field3D& vorticity_r,
    Field3D& vorticity_theta,
    Field3D& vorticity_z,
    Field3D& stretching_term,
    Field3D& tilting_term,
    Field3D& baroclinic_term)
{
    // Iterate over all grid points and compute vorticity diagnostics
    for (int i = 1; i < NR_ - 1; ++i) 
    {
        // Avoid division by zero by adding a small epsilon
        double r = i * dr_ + dynamics_constants::eps;

        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 1; k < NZ_ - 1; ++k) 
            {
                // Compute vorticity components
                double duth_dz = compute_dz(u_theta, i, j, k);
                double duz_dth = compute_dtheta(u_z, i, j, k);
                vorticity_r[i][j][k] = compute_vorticity_r(duz_dth, duth_dz);
                
                // Stretching term: ζ * ∂w/∂z
                double dw_dz = compute_dz(u_z, i, j, k);
                double zeta = vorticity_z[i][j][k];  // vertical vorticity
                stretching_term[i][j][k] = zeta * dw_dz;

                // Tilting term: ∂w/∂x * ∂v/∂z - ∂w/∂y * ∂u/∂z
                // In cylindrical coordinates: ∂w/∂r * ∂v/∂z - (1/r)∂w/∂θ * ∂u/∂z
                double dw_dr = compute_dr(u_z, i, j, k);
                double dv_dz = compute_dz(u_theta, i, j, k);
                double dw_dth = compute_dtheta(u_z, i, j, k);
                double du_dz = compute_dz(u_r, i, j, k);
                tilting_term[i][j][k] = dw_dr * dv_dz - (dw_dth / r) * du_dz;

                // Baroclinic generation: (1/ρ²) ∇ρ × ∇p
                double drho_dr = compute_dr(rho, i, j, k);
                double drho_dth = compute_dtheta(rho, i, j, k);
                double drho_dz = compute_dz(rho, i, j, k);

                double dp_dr = compute_dr(p, i, j, k);
                double dp_dth = compute_dtheta(p, i, j, k);
                double dp_dz = compute_dz(p, i, j, k);

                double rho_sq = rho[i][j][k] * rho[i][j][k];

                // Avoid division by zero by adding a small epsilon
                if (rho_sq > dynamics_constants::eps) 
                {
                    baroclinic_term[i][j][k] = (1.0 / rho_sq) *
                        (drho_dr * dp_dth - drho_dth * dp_dr); // z-component of cross product
                } 
                else 
                {
                    baroclinic_term[i][j][k] = 0.0;
                }
            }
        }
    }
}

/*This function computes the pressure diagnostics for the supercell scheme.
Takes in the u_r, u_theta, u_z, rho, theta, p_prime, dynamic_pressure, and buoyancy_pressure
and computes the pressure diagnostics for the supercell scheme.*/
void SupercellScheme::compute_pressure_diagnostics(
    const Field3D& u_r,
    const Field3D& u_theta,
    const Field3D& u_z,
    const Field3D& rho,
    const Field3D& theta,
    Field3D& p_prime,
    Field3D& dynamic_pressure,
    Field3D& buoyancy_pressure)
{
    // This would implement the Poisson equation solver for perturbation pressure
    // ∇²p' = -ρ₀ ∑ᵢⱼ ∂uᵢ/∂xⱼ ∂uⱼ/∂xᵢ + ∂(ρ₀ b)/∂z
    // For now, simplified perturbation pressure calculation

    // Iterate over all grid points and compute pressure diagnostics
    for (int i = 1; i < NR_ - 1; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 1; k < NZ_ - 1; ++k) 
            {
                // Simplified dynamic pressure (deformation/rotation component)
                double dur_dr = compute_dr(u_r, i, j, k);
                double dur_dth = compute_dtheta(u_r, i, j, k);
                double duth_dr = compute_dr(u_theta, i, j, k);
                double duth_dth = compute_dtheta(u_theta, i, j, k);
                double duz_dz = compute_dz(u_z, i, j, k);

                // Strain rate terms
                double deformation = dur_dr * dur_dr + (1.0/(i*dr_+dynamics_constants::eps)) *
                    (dur_dth * dur_dth + duth_dr * duth_dr + duth_dth * duth_dth) + duz_dz * duz_dz;

                dynamic_pressure[i][j][k] = -rho0_base[k] * deformation;

                // Buoyancy pressure (vertical buoyancy gradient)
                double theta_prime = theta[i][j][k] - theta0;
                double buoyancy = dynamics_constants::g * (theta_prime / theta0);
                buoyancy_pressure[i][j][k] = rho0_base[k] * buoyancy;

                // Total perturbation pressure
                p_prime[i][j][k] = dynamic_pressure[i][j][k] + buoyancy_pressure[i][j][k];
            }
        }
    }
}

// Helper function implementations

/*This function computes the derivative in the radial direction.
Takes in the field, the row index, the column index, and the level index
and computes the derivative in the radial direction.*/
double SupercellScheme::compute_dr(const Field3D& field, int i, int j, int k) const 
{
    return (field[i + 1][j][k] - field[i - 1][j][k]) / (2 * dr_);
}

/*This function computes the derivative in the azimuthal direction.
Takes in the field, the row index, the column index, and the level index
and computes the derivative in the azimuthal direction.*/
double SupercellScheme::compute_dtheta(const Field3D& field, int i, int j, int k) const 
{
    int j_prev = (j - 1 + NTH_) % NTH_;
    int j_next = (j + 1) % NTH_;
    return (field[i][j_next][k] - field[i][j_prev][k]) / (2 * dtheta_);
}

/*This function computes the derivative in the vertical direction.
Takes in the field, the row index, the column index, and the level index
and computes the derivative in the vertical direction.*/
double SupercellScheme::compute_dz(const Field3D& field, int i, int j, int k) const 
{
    return (field[i][j][k + 1] - field[i][j][k - 1]) / (2 * dz_);
}

/*This function computes the vorticity in the radial direction.
Takes in the dtheta_u_z, dz_u_theta, and computes the vorticity in the radial direction.*/
double SupercellScheme::compute_vorticity_r(double dtheta_u_z, double dz_u_theta) const 
{
    return dtheta_u_z - dz_u_theta;
}

/*This function computes the vorticity in the azimuthal direction.
Takes in the dz_u_r, dr_u_z, and computes the vorticity in the azimuthal direction.*/
double SupercellScheme::compute_vorticity_theta(double dz_u_r, double dr_u_z) const 
{
    return dz_u_r - dr_u_z;
}

/*This function computes the vorticity in the vertical direction.
Takes in the dr_u_theta, dtheta_u_r, and r
and computes the vorticity in the vertical direction.*/
double SupercellScheme::compute_vorticity_z(double dr_u_theta, double dtheta_u_r, double r) const 
{
    return (1.0 / r) * (dr_u_theta - dtheta_u_r) + dr_u_theta / r;
}

/*This function computes the buoyancy.
Takes in the theta_prime, rho, and rho0
and computes the buoyancy.*/
double SupercellScheme::compute_buoyancy(double theta_prime, double rho, double rho0) const 
{
    return dynamics_constants::g * (theta_prime / theta0) * (rho0 / rho);
}
