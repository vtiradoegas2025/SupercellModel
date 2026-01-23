# Source Code Structure

This directory contains the main implementation of SupercellModel's atmospheric simulation framework. The code is organized by physics modules following a factory pattern for extensibility.

## Core Files

### Main Simulation Components
- **`equations.cpp`** - Core atmospheric equations and field management
- **`dynamics.cpp`** - Dynamics module coordination and initialization
- **`numerics.cpp`** - Numerical methods coordination
- **`tornado_sim.cpp`** - Main simulation executable and configuration parsing
- **`gui.cpp`** - Optional SFML-based graphical interface

### Physics Module Coordinators
- **`boundary_layer.cpp`** - Planetary boundary layer coordination
- **`chaos.cpp`** - Stochastic perturbation coordination
- **`microphysics.cpp`** - Microphysics coordination
- **`radar.cpp`** - Radar forward operators coordination
- **`radiation.cpp`** - Radiation coordination
- **`terrain.cpp`** - Terrain/orography coordination
- **`turbulence.cpp`** - Sub-grid turbulence coordination

## Module Organization

Each physics module follows a consistent factory pattern:

```
module/
├── module.cpp          # Main coordination file
├── factory.cpp/.hpp    # Factory implementation for scheme selection
├── base/               # Common base classes and utilities
└── schemes/            # Individual physics schemes
    ├── scheme1/
    │   ├── scheme1.cpp
    │   └── scheme1.hpp
    └── scheme2/
        ├── scheme2.cpp
        └── scheme2.hpp
```

### boundary_layer/ - Planetary Boundary Layer
Parameterizes turbulent fluxes between surface and atmosphere.

**Schemes:**
- **`slab/`** - Simple slab model with constant flux layer
- **`ysu/`** - Yonsei University (YSU) nonlocal PBL scheme
- **`mynn/`** - Mellor-Yamada Nakanishi Niino (MYNN) TKE-based scheme

**Base components:**
- **`surface_fluxes.cpp/.hpp`** - Surface flux calculations (MOST similarity theory)

### chaos/ - Stochastic Perturbations
Implements controlled stochastic variability for ensemble forecasting.

**Schemes:**
- **`none/`** - No perturbations (deterministic)
- **`initial_conditions/`** - IC perturbations only
- **`boundary_layer/`** - BL tendency perturbations
- **`full_stochastic/`** - Complete SPPT implementation

**Base components:**
- **`random_generator.cpp/.hpp`** - Reproducible random number generation
- **`perturbation_field.cpp/.hpp`** - Field perturbation utilities
- **`correlation_filter.cpp/.hpp`** - Spatial correlation filtering

### dynamics/ - Large-Scale Dynamics
Handles compressible non-hydrostatic atmospheric dynamics.

**Schemes:**
- **`supercell/`** - Supercell-specific dynamics configuration
- **`tornado/`** - Tornado-scale dynamics with enhanced resolution

### microphysics/ - Cloud Microphysics
Parameterizes phase changes and precipitation processes.

**Schemes:**
- **`kessler/`** - Warm-rain Kessler scheme (autoconversion + accretion + evaporation)
- **`lin/`** - Lin et al. ice microphysics scheme
- **`thompson/`** - Thompson aerosol-aware mixed-phase scheme
- **`milbrandt/`** - Milbrandt-Yau double-moment scheme

**Base components:**
- **`thermodynamics.cpp/.hpp`** - Phase change thermodynamics and latent heating

### numerics/ - Numerical Methods
Core numerical algorithms for spatial discretization and time integration.

**Submodules:**
- **`advection/`** - Scalar transport schemes
  - **`tvd/`** - Total Variation Diminishing (monotonic)
  - **`weno5/`** - Weighted Essentially Non-Oscillatory (5th order)

- **`diffusion/`** - Momentum/heat diffusion
  - **`explicit/`** - Forward Euler diffusion
  - **`implicit/`** - Crank-Nicolson implicit diffusion

- **`time_stepping/`** - Time integration methods
  - **`rk3/`** - 3rd-order Runge-Kutta (Wicker & Skamarock)
  - **`rk4/`** - 4th-order Runge-Kutta

### radar/ - Radar Forward Operators
Forward operators for radar observables from model state.

**Schemes:**
- **`reflectivity/`** - Radar reflectivity (Z/Ze) from microphysics
- **`velocity/`** - Doppler radial velocity (Vr) from winds
- **`zdr/`** - Differential reflectivity (ZDR) from polarimetric variables

**Base components:**
- **`radar_base.cpp/.hpp`** - Common radar utilities and interfaces

### radiation/ - Atmospheric Radiation
Longwave and shortwave radiative transfer.

**Schemes:**
- **`simple_grey/`** - Gray atmosphere approximation with Newtonian cooling

**Base components:**
- **`radiative_transfer.cpp/.hpp`** - Radiative transfer utilities

### soundings/ - Atmospheric Soundings
Vertical profile initialization from observational data.

**Schemes:**
- **`sharpy/`** - SHARPY sounding format support

**Base components:**
- **`soundings_base.cpp/.hpp`** - Sounding data structures and interpolation

### terrain/ - Terrain/Orography
Topographic effects on atmospheric flow.

**Schemes:**
- **`none/`** - Flat terrain
- **`bell/`** - Axisymmetric bell-shaped mountain
- **`schar/`** - Schar-type ridge topography

**Base components:**
- **`topography.cpp/.hpp`** - Terrain generation and coordinate transformation

### turbulence/ - Sub-Grid Turbulence
Parameterizes unresolved turbulent motions.

**Schemes:**
- **`smagorinsky/`** - Constant Smagorinsky-Lilly eddy viscosity
- **`tke/`** - Prognostic turbulence kinetic energy scheme

**Base components:**
- **`eddy_viscosity.cpp/.hpp`** - Eddy viscosity calculations

## Development Notes

### Code Style
- **C++17 standard** with modern features
- **RAII principles** for resource management
- **Factory pattern** for physics scheme selection
- **Header-only base classes** where appropriate
- **Comprehensive error handling** with descriptive messages

### Build System
- **Makefile-based** build with automatic dependency detection
- **Modular compilation** - individual modules can be excluded
- **SFML integration** for optional GUI (controlled by `GUI=1` make variable)

### Testing
- **Integration tests** in `tests/test_integration.sh`
- **Sounding tests** in `soundings/test_soundings.cpp`
- **Manual verification** through visualization pipeline

### Extension Points
- **New physics schemes**: Add to appropriate `schemes/` subdirectory
- **New modules**: Follow existing factory pattern structure
- **Configuration**: Extend YAML parsing in `tornado_sim.cpp`

### Performance Considerations
- **Cylindrical coordinates** for efficient azimuthal symmetry
- **Factory initialization** done once at startup
- **Memory pools** for frequently allocated objects
- **SIMD-friendly data layouts** where beneficial

This modular structure enables easy comparison of physics schemes, simplified testing, and straightforward extension of model capabilities.
