#pragma once
#include <vector>
#include <memory>
#include <string>
#include "numerics_base.hpp"

/*This header file contains the base classes and structures for the diffusion module.
The diffusion module is responsible for the diffusion of the scalar fields in the simulation.
The diffusion scheme is chosen by the user in the configuration file.
This module is used to compute the diffusion of the simulation.*/


// Diffusion module configuration
struct DiffusionConfig 
{
    std::string scheme_id = "explicit";
    std::string operator_type = "laplacian";     // "laplacian" or "anisotropic"
    std::string apply_to = "scalars";            // "scalars", "momentum", "all"
    double K_h = 0.0;                           // horizontal diffusivity [m²/s]
    double K_v = 0.0;                           // vertical diffusivity [m²/s]
    std::string implicit_dim = "none";           // "none", "vertical_only", "full_3d"
    double dt_diffusion = 1.0;                   // diffusion timestep [s]
    bool use_variable_K = false;                 // use variable diffusivity fields
};

// State view for diffusion computations
struct DiffusionStateView {
    // Fields to diffuse
    const std::vector<std::vector<std::vector<double>>>* u = nullptr;      // zonal momentum
    const std::vector<std::vector<std::vector<double>>>* v = nullptr;      // meridional momentum
    const std::vector<std::vector<std::vector<double>>>* w = nullptr;      // vertical momentum
    const std::vector<std::vector<std::vector<double>>>* theta = nullptr;  // potential temperature
    const std::vector<std::vector<std::vector<double>>>* qv = nullptr;     // water vapor
    const std::vector<std::vector<std::vector<double>>>* rho = nullptr;    // density

    // Diffusivity fields (if use_variable_K = true)
    const std::vector<std::vector<std::vector<double>>>* K_momentum = nullptr;  // momentum diffusivity
    const std::vector<std::vector<std::vector<double>>>* K_scalar = nullptr;    // scalar diffusivity

    // Grid information
    const GridMetrics* grid = nullptr;
};

// Output tendencies from diffusion
struct DiffusionTendencies {
    std::vector<std::vector<std::vector<double>>> dudt_diff;     // u diffusion tendency
    std::vector<std::vector<std::vector<double>>> dvdt_diff;     // v diffusion tendency
    std::vector<std::vector<std::vector<double>>> dwdt_diff;     // w diffusion tendency
    std::vector<std::vector<std::vector<double>>> dthetadt_diff; // theta diffusion tendency
    std::vector<std::vector<std::vector<double>>> dqvdt_diff;    // qv diffusion tendency
};

// Diagnostics for diffusion schemes
struct DiffusionDiagnostics 
{
    double max_diffusion_number = 0.0;  // maximum diffusion number for stability
    std::vector<std::vector<std::vector<double>>> K_effective;  // effective diffusivity field
};

// Abstract base class for diffusion schemes
class DiffusionSchemeBase : public NumericalSchemeBase 
{
public:
    virtual void initialize(const DiffusionConfig& cfg) = 0;

    virtual void compute_diffusion_tendencies(
        const DiffusionConfig& cfg,
        const DiffusionStateView& state,
        DiffusionTendencies& tendencies,
        DiffusionDiagnostics* diag_opt = nullptr) = 0;

    virtual double check_stability(
        const DiffusionConfig& cfg,
        const DiffusionStateView& state) = 0;
};

// Factory function type
using DiffusionSchemeFactory = std::unique_ptr<DiffusionSchemeBase> (*)(const DiffusionConfig&);

// Global diffusion scheme instance and configuration declared in simulation.hpp

// Factory and initialization functions
std::unique_ptr<DiffusionSchemeBase> create_diffusion_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_diffusion_schemes();
void initialize_diffusion(const std::string& scheme_name,
                         const DiffusionConfig& cfg = DiffusionConfig{});
