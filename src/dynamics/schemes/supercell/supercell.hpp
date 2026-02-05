#pragma once
#include "dynamics_base.hpp"
#include <vector>


/*This class is the supercell scheme*/
class SupercellScheme : public DynamicsScheme
{
public:
    SupercellScheme();

    /*This function computes the momentum tendencies for the supercell scheme.
    Takes in the u_r, u_theta, u_z, rho, p, theta, dt, du_r_dt, du_theta_dt, du_z_dt, drho_dt, and dp_dt
    and computes the momentum tendencies for the supercell scheme.*/
    void compute_momentum_tendencies(
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
        Field3D& dp_dt
    ) override;

    /*This function computes the vorticity diagnostics for the supercell scheme.
    Takes in the u_r, u_theta, u_z, rho, p, vorticity_r, vorticity_theta, vorticity_z, stretching_term, tilting_term, and baroclinic_term
    and computes the vorticity diagnostics for the supercell scheme.*/
    void compute_vorticity_diagnostics(
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
        Field3D& baroclinic_term
    ) override;

    /*This function computes the pressure diagnostics for the supercell scheme.
    Takes in the u_r, u_theta, u_z, rho, theta, p_prime, dynamic_pressure, and buoyancy_pressure
    and computes the pressure diagnostics for the supercell scheme.*/
    void compute_pressure_diagnostics(
        const Field3D& u_r,
        const Field3D& u_theta,
        const Field3D& u_z,
        const Field3D& rho,
        const Field3D& theta,
        Field3D& p_prime,
        Field3D& dynamic_pressure,
        Field3D& buoyancy_pressure
    ) override;

    std::string get_scheme_name() const override { return "supercell"; }
    std::string get_coordinate_system() const override { return "cylindrical"; }
    int get_num_prognostic_vars() const override { return 5; } // u_r, u_theta, u_z, rho, p

private:
    // Helper functions for computing derivatives in cylindrical coordinates
    double compute_dr(const Field3D& field, int i, int j, int k) const;
    double compute_dtheta(const Field3D& field, int i, int j, int k) const;
    double compute_dz(const Field3D& field, int i, int j, int k) const;

    // Vorticity component calculations
    double compute_vorticity_r(double dtheta_u_z, double dz_u_theta) const;
    double compute_vorticity_theta(double dz_u_r, double dr_u_z) const;
    double compute_vorticity_z(double dr_u_theta, double dtheta_u_r, double r) const;

    // Buoyancy calculation
    double compute_buoyancy(double theta_prime, double rho, double rho0) const;

    // Grid parameters (will be set from global simulation parameters)
    int NR_, NTH_, NZ_;
    double dr_, dtheta_, dz_;
};
