#pragma once
#include "dynamics_base.hpp"
#include <vector>


/*This header file contains the declaration of the TornadoScheme class.
This class implements the tornado scheme.
It is a subclass of the DynamicsScheme class.
It implements the compute_momentum_tendencies method.
It implements the compute_angular_momentum method.
It implements the compute_vorticity_diagnostics method.
*/
class TornadoScheme : public DynamicsScheme
{
public:
    TornadoScheme();

    // Core momentum tendencies (axisymmetric dynamics)
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

    // Angular momentum diagnostics (key for tornado spin-up)
    void compute_angular_momentum(
        const std::vector<std::vector<std::vector<float>>>& u_r,
        const std::vector<std::vector<std::vector<float>>>& u_theta,
        std::vector<std::vector<std::vector<float>>>& angular_momentum,
        std::vector<std::vector<std::vector<float>>>& angular_momentum_tendency
    ) override;

    // Axisymmetric vorticity diagnostics
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

    std::string get_scheme_name() const override { return "tornado"; }
    std::string get_coordinate_system() const override { return "cylindrical"; }
    int get_num_prognostic_vars() const override { return 5; } // u_r, u_theta, u_z, rho, p

private:
    // Helper functions for axisymmetric calculations
    double compute_dr(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const;
    double compute_dz(const std::vector<std::vector<std::vector<float>>>& field, int i, int j, int k) const;

    // Radial mass flux diagnostic (for vertical motion)
    double compute_radial_mass_flux(const std::vector<std::vector<std::vector<float>>>& u_r,
                                   const std::vector<std::vector<std::vector<float>>>& rho,
                                   int i, int k) const;

    // Grid parameters
    int NR_, NTH_, NZ_;
    double dr_, dtheta_, dz_;
};
