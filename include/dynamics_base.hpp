#pragma once
#include <vector>
#include <memory>
#include <string>

/*This header file contains the base classes and structures for the dynamics module.
The dynamics module is responsible for the dynamics of the simulation.
The dynamics scheme is chosen by the user in the configuration file.
This module is used to compute the dynamics of the simulation.*/



class DynamicsScheme
{
public:
    virtual ~DynamicsScheme() = default;

    // Core interface - compute tendencies for momentum equations
    virtual void compute_momentum_tendencies(
        // Input state (cylindrical coordinates: r, θ, z)
        const std::vector<std::vector<std::vector<float>>>& u_r,        // radial wind [r][θ][z]
        const std::vector<std::vector<std::vector<float>>>& u_theta,    // azimuthal wind [r][θ][z]
        const std::vector<std::vector<std::vector<float>>>& u_z,        // vertical wind [r][θ][z]
        const std::vector<std::vector<std::vector<float>>>& rho,        // density [r][θ][z]
        const std::vector<std::vector<std::vector<float>>>& p,          // pressure [r][θ][z]
        const std::vector<std::vector<std::vector<float>>>& theta,      // potential temperature [r][θ][z]
        double dt,                                                      // timestep
        // Output tendencies (will be added to existing fields)
        std::vector<std::vector<std::vector<float>>>& du_r_dt,          // radial momentum tendency
        std::vector<std::vector<std::vector<float>>>& du_theta_dt,      // azimuthal momentum tendency
        std::vector<std::vector<std::vector<float>>>& du_z_dt,          // vertical momentum tendency
        std::vector<std::vector<std::vector<float>>>& drho_dt,          // density tendency
        std::vector<std::vector<std::vector<float>>>& dp_dt             // pressure tendency
    ) = 0;

    // Vorticity diagnostics (key for supercell/tornado development)
    virtual void compute_vorticity_diagnostics(
        const std::vector<std::vector<std::vector<float>>>& u_r,
        const std::vector<std::vector<std::vector<float>>>& u_theta,
        const std::vector<std::vector<std::vector<float>>>& u_z,
        const std::vector<std::vector<std::vector<float>>>& rho,
        const std::vector<std::vector<std::vector<float>>>& p,
        // Output vorticity components and tendencies
        std::vector<std::vector<std::vector<float>>>& vorticity_r,      // radial vorticity
        std::vector<std::vector<std::vector<float>>>& vorticity_theta,  // azimuthal vorticity
        std::vector<std::vector<std::vector<float>>>& vorticity_z,      // vertical vorticity (ζ)
        std::vector<std::vector<std::vector<float>>>& stretching_term,  // ζ * ∂w/∂z
        std::vector<std::vector<std::vector<float>>>& tilting_term,     // tilting of horizontal vorticity
        std::vector<std::vector<std::vector<float>>>& baroclinic_term   // baroclinic generation
    ) {}

    // Angular momentum diagnostics (for tornado dynamics)
    virtual void compute_angular_momentum(
        const std::vector<std::vector<std::vector<float>>>& u_r,
        const std::vector<std::vector<std::vector<float>>>& u_theta,
        std::vector<std::vector<std::vector<float>>>& angular_momentum,  // M = r * V
        std::vector<std::vector<std::vector<float>>>& angular_momentum_tendency
    ) {}

    // Pressure diagnostics (perturbation pressure solver for supercells)
    virtual void compute_pressure_diagnostics(
        const std::vector<std::vector<std::vector<float>>>& u_r,
        const std::vector<std::vector<std::vector<float>>>& u_theta,
        const std::vector<std::vector<std::vector<float>>>& u_z,
        const std::vector<std::vector<std::vector<float>>>& rho,
        const std::vector<std::vector<std::vector<float>>>& theta,
        std::vector<std::vector<std::vector<float>>>& p_prime,          // perturbation pressure
        std::vector<std::vector<std::vector<float>>>& dynamic_pressure,  // deformation/rotation component
        std::vector<std::vector<std::vector<float>>>& buoyancy_pressure  // buoyancy component
    ) {}

    // Get scheme name for diagnostics
    virtual std::string get_scheme_name() const = 0;

    // Get coordinate system used (cartesian/cylindrical)
    virtual std::string get_coordinate_system() const = 0;

    // Get number of prognostic variables
    virtual int get_num_prognostic_vars() const = 0;
};

// Factory function declaration
std::unique_ptr<DynamicsScheme> create_dynamics_scheme(const std::string& scheme_name);

// Common physical constants used by dynamics schemes
namespace dynamics_constants
{
    constexpr double g = 9.81;           // gravity (m/s²)
    constexpr double f_coriolis = 1e-4;  // Coriolis parameter (s⁻¹) - often neglected
    constexpr double R_d = 287.0;        // dry air gas constant (J/kg·K)
    constexpr double cp = 1004.0;        // specific heat at constant pressure (J/kg·K)
    constexpr double p0 = 100000.0;      // reference pressure (Pa)
    constexpr double theta0 = 300.0;     // base potential temperature (K)
    constexpr double gamma = cp / (cp - R_d); // adiabatic index
    constexpr double L_v = 2.5e6;        // latent heat of vaporization (J/kg)

    // Numerical constants
    constexpr double eps = 1e-8;         // small number for avoiding division by zero
    constexpr double nu_visc = 0.01;     // numerical viscosity coefficient
    constexpr double sponge_strength = 0.01; // Rayleigh sponge damping
}
