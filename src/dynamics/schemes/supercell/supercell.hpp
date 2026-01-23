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
        std::vector<std::vector<std::vector<float>>>& dp_dt
    ) override;

    /*This function computes the vorticity diagnostics for the supercell scheme.
    Takes in the u_r, u_theta, u_z, rho, p, vorticity_r, vorticity_theta, vorticity_z, stretching_term, tilting_term, and baroclinic_term
    and computes the vorticity diagnostics for the supercell scheme.*/
    void compute_vorticity_diagnostics(
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
        std::vector<std::vector<std::vector<float>>>& baroclinic_term
    ) override;

    /*This function computes the pressure diagnostics for the supercell scheme.
    Takes in the u_r, u_theta, u_z, rho, theta, p_prime, dynamic_pressure, and buoyancy_pressure
    and computes the pressure diagnostics for the supercell scheme.*/
    void compute_pressure_diagnostics(
        const std::vector<std::vector<std::vector<float>>>& u_r,
        const std::vector<std::vector<std::vector<float>>>& u_theta,
        const std::vector<std::vector<std::vector<float>>>& u_z,
        const std::vector<std::vector<std::vector<float>>>& rho,
        const std::vector<std::vector<std::vector<float>>>& theta,
        std::vector<std::vector<std::vector<float>>>& p_prime,
        std::vector<std::vector<std::vector<float>>>& dynamic_pressure,
        std::vector<std::vector<std::vector<float>>>& buoyancy_pressure
    ) override;

    std::string get_scheme_name() const override { return "supercell"; }
    std::string get_coordinate_system() const override { return "cylindrical"; }
    int get_num_prognostic_vars() const override { return 5; } // u_r, u_theta, u_z, rho, p

private:
    // Helper functions for computing derivatives in cylindrical coordinates
    double compute_dr(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const;
    double compute_dtheta(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const;
    double compute_dz(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const;

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
