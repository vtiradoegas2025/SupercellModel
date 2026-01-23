#pragma once
#include <vector>
#include <memory>
#include <string>
#include "numerics_base.hpp"


/*This header file contains the base classes and structures for the time stepping module.
The time stepping module is responsible for the time stepping of the simulation.
The time stepping scheme is chosen by the user in the configuration file.
This module is used to time step the simulation and ensure the stability of the simulation.*/

// Time stepping module configuration
struct TimeSteppingConfig 
{
    std::string scheme_id = "rk3";
    bool split_acoustic = false;              // split-explicit compressible
    int n_acoustic_substeps = 1;              // acoustic substeps if split
    std::string physics_splitting = "additive"; // "additive" or "strang"
    double cfl_safety = numerics_constants::cfl_target;
    bool adaptive_dt = false;                 // adaptive timestepping
    double dt_min = 1e-6;                     // minimum timestep [s]
    double dt_max = 100.0;                    // maximum timestep [s]
};

// State container for time stepping
struct TimeSteppingState {
    std::vector<NumericalState> fields;       // all prognostic fields
    double time = 0.0;                        // current time [s]
    double dt = 1.0;                          // current timestep [s]
};

// Diagnostics for time stepping
struct TimeSteppingDiagnostics {
    int n_steps = 0;                          // number of steps taken
    double max_cfl = 0.0;                     // maximum CFL number
    double cpu_time_per_step = 0.0;           // performance metric
    std::vector<double> dt_history;            // timestep history
    bool converged = true;                    // convergence flag
};

// Abstract base class for time stepping schemes
class TimeSteppingSchemeBase : public NumericalSchemeBase {
public:
    virtual void initialize(const TimeSteppingConfig& cfg, RHSFunction rhs_func) = 0;

    virtual void step(
        const TimeSteppingConfig& cfg,
        TimeSteppingState& state,
        TimeSteppingDiagnostics* diag_opt = nullptr) = 0;

    virtual double suggest_dt(
        const TimeSteppingConfig& cfg,
        const TimeSteppingState& state) = 0;
};

// Factory function type
using TimeSteppingSchemeFactory = std::unique_ptr<TimeSteppingSchemeBase> (*)(const TimeSteppingConfig&);

// Global time stepping scheme instance and configuration declared in simulation.hpp

// Factory and initialization functions
std::unique_ptr<TimeSteppingSchemeBase> create_time_stepping_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_time_stepping_schemes();
void initialize_time_stepping(const std::string& scheme_name,
                             const TimeSteppingConfig& cfg = TimeSteppingConfig{},
                             RHSFunction rhs_func = nullptr);
