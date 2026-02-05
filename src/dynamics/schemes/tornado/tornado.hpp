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

    // Angular momentum diagnostics (key for tornado spin-up)
    void compute_angular_momentum(
        const Field3D& u_r,
        const Field3D& u_theta,
        Field3D& angular_momentum,
        Field3D& angular_momentum_tendency
    ) override;

    // Axisymmetric vorticity diagnostics
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

    std::string get_scheme_name() const override { return "tornado"; }
    std::string get_coordinate_system() const override { return "cylindrical"; }
    int get_num_prognostic_vars() const override { return 5; } // u_r, u_theta, u_z, rho, p

private:
    // Helper functions for axisymmetric calculations
    double compute_dr(const Field3D& field, int i, int j, int k) const;
    double compute_dz(const Field3D& field, int i, int j, int k) const;

    // Radial mass flux diagnostic (for vertical motion)
    double compute_radial_mass_flux(const Field3D& u_r,
                                   const Field3D& rho,
                                   int i, int k) const;

    // Grid parameters
    int NR_, NTH_, NZ_;
    double dr_, dtheta_, dz_;
};
