#pragma once
#include <vector>
#include <memory>
#include <string>
#include "numerics_base.hpp"

/*This header file contains the base classes and structures for the advection module.
The advection module is responsible for advecting the scalar fields in the simulation.
The advection scheme is chosen by the user in the configuration file.
This module is used to compute the advection of the simulation.*/

/*Advection module configuration*/
struct AdvectionConfig 
{
    std::string scheme_id = "tvd";
    std::string formulation = "flux_form";  // "flux_form" or "advective_form"
    std::string split_mode = "unsplit";     // "unsplit", "strang_split", "directional_split"
    std::string limiter_id = "mc";          // TVD limiter: "minmod", "vanleer", "superbee", "mc", "universal"
    bool positivity = false;                // enforce q >= 0 for tracers
    double cfl_target = numerics_constants::cfl_target;
    std::string reconstruct_on = "q";       // "q" or "rhoq"
    int weno_order = 5;                     // WENO order (5 for WENO5)
};

/*State view for advection computations*/
struct AdvectionStateView 
{
    // Velocity fields (cell-centered)
    const std::vector<std::vector<std::vector<double>>>* u = nullptr;    // zonal wind [m/s]
    const std::vector<std::vector<std::vector<double>>>* v = nullptr;    // meridional wind [m/s]
    const std::vector<std::vector<std::vector<double>>>* w = nullptr;    // vertical wind [m/s]

    // Tracer fields to advect
    const std::vector<std::vector<std::vector<double>>>* q = nullptr;    // scalar field(s)
    const std::vector<std::vector<std::vector<double>>>* rho = nullptr;  // density (optional)

    // Grid information
    const GridMetrics* grid = nullptr;
};

// Output tendencies from advection
struct AdvectionTendencies 
{
    std::vector<std::vector<std::vector<double>>> dqdt_adv;  // advection tendency [units/s]
};

// Diagnostics for advection schemes
struct AdvectionDiagnostics 
{
    double max_cfl_x = 0.0;     // maximum CFL in x-direction
    double max_cfl_y = 0.0;     // maximum CFL in y-direction
    double max_cfl_z = 0.0;     // maximum CFL in z-direction
    double suggested_dt = 0.0;  // suggested timestep for stability
    std::vector<std::vector<std::vector<double>>> fluxes_x;  // x-fluxes (optional)
    std::vector<std::vector<std::vector<double>>> fluxes_y;  // y-fluxes (optional)
    std::vector<std::vector<std::vector<double>>> fluxes_z;  // z-fluxes (optional)
};

/*Abstract base class for advection schemes*/
class AdvectionSchemeBase : public NumericalSchemeBase 
{
public:
    virtual void initialize(const AdvectionConfig& cfg) = 0;

    //Compute the flux divergence of the scalar fields
    virtual void compute_flux_divergence(
        const AdvectionConfig& cfg,
        const AdvectionStateView& state,
        AdvectionTendencies& tendencies,
        AdvectionDiagnostics* diag_opt = nullptr) = 0;

    virtual double suggest_dt(
        const AdvectionConfig& cfg,
        const AdvectionStateView& state) = 0;
};

// Factory function type
using AdvectionSchemeFactory = std::unique_ptr<AdvectionSchemeBase> (*)(const AdvectionConfig&);

// Global advection scheme instance and configuration declared in simulation.hpp

// Factory and initialization functions
std::unique_ptr<AdvectionSchemeBase> create_advection_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_advection_schemes();
void initialize_advection(const std::string& scheme_name,
                         const AdvectionConfig& cfg = AdvectionConfig{});
