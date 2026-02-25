# Header Files

This directory contains the header files defining interfaces, base classes, and type definitions for SupercellModel. Headers are organized by physics module and provide the public API for the simulation framework.

## Core Headers

### Main Simulation Interface
- **`simulation.hpp`** - Primary simulation class and configuration structures
  - Defines `SimulationState`, `WindProfile`, grid metrics
  - Contains main initialization and stepping functions
  - Entry point for external coupling

### Runtime/Validation Headers
- **`runtime_config.hpp`** - Runtime configuration globals and parser-facing declarations
- **`field_contract.hpp`** - CM1-style export/validation contract metadata
- **`field_validation.hpp`** - Field-level guard/validation result structures and checks
- **`physical_constants.hpp`** - Shared constants consumed across physics modules

## Physics Module Base Classes

Each physics module provides a base class defining the interface for all schemes within that module. These headers contain:

### advection_base.hpp - Scalar Transport
- **`AdvectionScheme`** - Base class for advection algorithms
- **TVD/MC and WENO implementations**
- **Conservative transport interfaces**

### boundary_layer_base.hpp - Planetary Boundary Layer
- **`BoundaryLayerScheme`** - PBL parameterization interface
- **Surface flux calculations (MOST similarity theory)**
- **Non-local mixing algorithms**
- **See `src/boundary_layer/README.md` for implementation details**

### chaos_base.hpp - Stochastic Perturbations
- **`ChaosScheme`** - Stochastic parameterization interface
- **Perturbation field definitions**
- **SPPT (Stochastically Perturbed Parameterization Tendencies)**
- **See `src/chaos/README.md` for implementation details**

### diffusion_base.hpp - Momentum/Heat Diffusion
- **`DiffusionScheme`** - Diffusion algorithm interface
- **Explicit and implicit Laplacian methods**
- **Hyper-viscosity options**
- **See `src/numerics/README.md` for implementation details**

### dynamics_base.hpp - Large-Scale Dynamics
- **`DynamicsScheme`** - Atmospheric dynamics interface
- **Compressible Euler equations**
- **Pressure gradient and buoyancy terms**
- **See `src/dynamics/README.md` for implementation details**

### microphysics_base.hpp - Cloud Microphysics
- **`MicrophysicsScheme`** - Microphysics parameterization interface
- **Phase change rates and latent heating**
- **Precipitation fallout algorithms**
- **See `src/microphysics/README.md` for implementation details**

### numerics_base.hpp - Numerical Methods Coordination
- **`NumericsManager`** - Coordinator for advection/diffusion/time-stepping
- **CFL condition monitoring**
- **Stability diagnostics**
- **See `src/numerics/README.md` for implementation details**

### radar_base.hpp - Radar Forward Operators
- **`RadarSchemeBase`** - Base class for radar observation operators
- **`RadarConfig`** - Configuration structure for radar settings
- **`RadarStateView`** - Read-only view of model state for radar calculations
- **`RadarOut`** - Output structure for radar observables
- **Reflectivity operators** (Z/Ze from hydrometeor mixing ratios)
- **Doppler velocity operators** (V_r from wind field projections)
- **Polarimetric operators** (Z_DR with ZH/ZV components)
- **See `src/radar/README.md` for implementation details**

### radiation_base.hpp - Atmospheric Radiation
- **`RadiationScheme`** - Radiative transfer interface
- **Longwave/shortwave cooling/heating**
- **Gray atmosphere approximations**
- **See `src/radiation/README.md` for implementation details**

### soundings_base.hpp - Atmospheric Soundings
- **`SoundingScheme`** - Vertical profile interface
- **Observational data ingestion**
- **Profile interpolation algorithms**
- **See `src/soundings/README.md` for implementation details**

### soundings.hpp - Sounding Data Structures
- **Sounding data containers**
- **Profile manipulation utilities**
- **Quality control interfaces**
- **See `src/soundings/README.md` for implementation details**

### terrain_base.hpp - Terrain/Orography
- **`TerrainScheme`** - Terrain parameterization interface
- **Topographic forcing**
- **Coordinate transformation**
- **See `src/terrain/README.md` for implementation details**

### time_stepping_base.hpp - Time Integration
- **`TimeSteppingScheme`** - Time integration interface
- **Runge-Kutta implementations**
- **Splitting methods (HEVI)**
- **See `src/numerics/README.md` for implementation details**

### turbulence_base.hpp - Sub-Grid Turbulence
- **`TurbulenceScheme`** - Turbulence parameterization interface
- **Eddy viscosity models**
- **TKE prognostic schemes**
- **See `src/turbulence/README.md` for implementation details**

## Design Principles

### Interface Consistency
All base classes follow consistent patterns:
- **Factory registration** via `register_*_scheme()` functions
- **Initialization** through `initialize()` methods
- **Application** via `apply_*()` or `compute_*()` methods
- **Configuration** through parameter structures

### Memory Management
- **RAII principles** throughout
- **Smart pointer usage** where appropriate
- **No raw pointers** in public interfaces
- **Exception safety** guaranteed

### Type Safety
- **Strong typing** with custom structs/enums
- **Template metaprogramming** for compile-time optimization
- **Runtime type checking** with assertions
- **Unit-safe quantities** where beneficial

### Extensibility
- **Plugin architecture** via factory registration
- **Forward-compatible interfaces**
- **Versioned base classes**
- **Optional feature flags**

## Usage Patterns

### Implementing New Schemes
```cpp
#include "microphysics_base.hpp"

class MyMicrophysics : public MicrophysicsScheme {
public:
    void initialize(const MicrophysicsConfig& config) override {
        // Implementation
    }

    void compute_tendencies(/* parameters */) override {
        // Implementation
    }
};

// Register with factory
REGISTER_MICROPHYSICS_SCHEME(MyMicrophysics, "my_scheme");
```

### Using Physics Modules
```cpp
#include "simulation.hpp"
#include "microphysics_base.hpp"

// Initialize simulation
SimulationState state;
initialize_simulation(state, config);

// Get microphysics scheme
auto microphysics = create_microphysics_scheme("kessler");
microphysics->initialize(microphysics_config);

// Apply in time loop
microphysics->compute_tendencies(state, tendencies, dt);
```

## Build System Integration

Headers are included via the main Makefile:
- **Automatic dependency generation**
- **Precompiled headers** support
- **Include path management**
- **Conditional compilation** flags

## Testing

Headers include extensive contracts:
- **Preconditions** checked via assertions
- **Postconditions** validated
- **Invariant maintenance**
- **Error handling** with descriptive messages

## Documentation

Each header contains:
- **Doxygen-compatible comments**
- **Mathematical formulations**
- **Algorithm references**
- **Usage examples**
- **Performance characteristics**

This header organization provides a clean separation between interface and implementation, enabling modular development and testing of atmospheric physics parameterizations.
