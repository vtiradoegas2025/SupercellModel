#pragma once
#include <vector>
#include <memory>
#include <string>
#include "turbulence_base.hpp"
#include <functional>

/*This header file contains the base classes and structures for the numerics module.
The numerics module is responsible for the numerics of the simulation.
The numerics scheme is chosen by the user in the configuration file.
This module is used to compute the numerics of the simulation.*/


// Numerics module constants
namespace numerics_constants 
{
    inline constexpr double epsilon = 1e-12;        // small number for divisions
    inline constexpr double weno_epsilon = 1e-6;    // WENO epsilon for smoothness
    inline constexpr double cfl_target = 0.8;       // target CFL number
}

// GridMetrics is defined in turbulence_base.hpp

// Generic state vector for numerical operations
struct NumericalState 
{
    std::vector<std::vector<std::vector<double>>> data;  // 3D field data
    std::string name;  // field name for diagnostics
};

// Generic tendencies for numerical operations
struct NumericalTendencies 
{
    std::vector<std::vector<std::vector<double>>> tendencies;  // 3D tendency data
    std::string name;  // field name
};

// RHS function type for time stepping
using RHSFunction = std::function<void(const std::vector<NumericalState>&, double,
                                       std::vector<NumericalTendencies>&)>;

// Base interface for all numerical schemes
class NumericalSchemeBase 
{
public:
    virtual ~NumericalSchemeBase() = default;
    virtual std::string name() const = 0;
    virtual void initialize() = 0;
};

// Global numerics configuration (set from YAML)
extern GridMetrics global_grid_metrics;
