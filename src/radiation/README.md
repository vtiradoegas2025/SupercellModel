# Atmospheric Radiation Module

This module implements radiative transfer parameterizations for longwave and shortwave radiation, providing heating/cooling rates that drive atmospheric temperature changes.

## Overview

Atmospheric radiation affects:
- **Atmospheric stability**: Radiative heating/cooling profiles
- **Boundary layer development**: Surface and nocturnal cooling
- **Cloud formation**: Radiative-convective equilibrium
- **Diurnal cycles**: Day-night temperature variations

## Architecture

```
src/radiation/
├── radiation.cpp            # Module coordinator
├── factory.cpp/.hpp         # Radiation scheme factory
├── base/
│   └── radiative_transfer.cpp/.hpp # Radiative transfer utilities
└── schemes/                 # Radiation implementations
    └── simple_grey/         # Grey atmosphere approximation
```

## Radiation Physics

### Radiative Transfer Equation

The fundamental equation for radiative intensity I_ν:

```
dI_ν/ds = -κ_ν ρ I_ν + κ_ν ρ B_ν(T) + S_ν
```

Where:
- I_ν: Specific intensity
- κ_ν: Absorption coefficient
- ρ: Density
- B_ν(T): Planck function
- S_ν: Source terms

### Heating Rate

Net radiative heating/cooling rate:

```
Q_radiation = (1/ρ c_p) ∇ · F
```

Where F is the net radiative flux.

## Scheme Implementation

### Simple Grey Atmosphere

**Assumptions:**
- Single broadband absorption coefficient
- No wavelength dependence
- Newtonian relaxation to reference profile

**Longwave Radiation:**
```cpp
// LW cooling rate
dT/dt = - (σ T^4) / (ρ c_p H_p)
```

Where:
- σ: Stefan-Boltzmann constant
- H_p: Pressure scale height
- ρ: Density
- c_p: Specific heat

**Shortwave Radiation:**
```cpp
// SW heating (simplified)
dT/dt = S_0 (1 - α) exp(-τ_z / μ) / (ρ c_p H_p)
```

Where:
- S_0: Solar constant
- α: Surface albedo
- τ_z: Optical depth
- μ: Cosine of zenith angle

### RRTMG (Rapid Radiative Transfer Model for GCMs)

**Features (when implemented):**
- 16 spectral bands (LW), 14 bands (SW)
- Correlated-k distribution method
- Multiple scattering calculations
- Cloud optical properties

## Implementation Details

### Radiative Transfer Solver

**Two-stream approximation:**
```cpp
// Upward and downward fluxes
dF↑/dz = -κ ρ (F↑ - π B)
dF↓/dz = κ ρ (F↓ - π B)
```

**Solution methods:**
- **Adding method**: Exact analytical solution
- **Discrete ordinates**: Numerical integration
- **Eddington approximation**: Simplified two-stream

### Time Integration

**Coupling with dynamics:**
```cpp
// Radiative tendency
dθ/dt = Q_radiation / (c_p π_ref^{R_d/c_p})
```

**Time stepping:**
```cpp
// Explicit coupling
θ^{n+1} = θ^n + Δt * Q_radiation(θ^n, q^n)

// Implicit coupling (stable)
(1 - Δt ∂Q/∂θ) θ^{n+1} = θ^n + Δt Q_radiation(q^n)
```

### Cloud Effects

**Cloud optical properties:**
```cpp
// Liquid water optical depth
τ_cloud = ∫ κ_lwc q_l dz
```

Where κ_lwc depends on droplet size distribution.

## Configuration

### Basic Setup
```yaml
radiation:
  scheme: "simple_grey"     # simple_grey, rrtmg
  dt_radiation: 300.0       # Radiation time step (s)
  do_lw: true               # Longwave radiation
  do_sw: true               # Shortwave radiation
```

### Grey Atmosphere Parameters
```yaml
radiation:
  # LW parameters
  tau_lw_ref: 0.1           # Reference LW optical depth
  lw_relaxation_time: 3600.0 # Relaxation time (s)

  # SW parameters
  solar_constant: 1366.0    # Solar constant (W/m²)
  surface_albedo: 0.2       # Surface reflectivity
  sw_optical_depth: 0.05    # SW optical depth
```

### Advanced Options
```yaml
radiation:
  # Cloud effects
  cloud_lw_overlap: "random"  # Cloud overlap assumption
  cloud_sw_scattering: true   # Include scattering

  # Numerical options
  n_streams: 2               # Number of streams (2, 4, 8)
  solver_tolerance: 1e-6     # Iterative solver tolerance
```

## Validation and Testing

### Radiative Equilibrium Tests
```bash
# Test radiative-convective equilibrium
python tests/test_radiation.py --scheme simple_grey --equilibrium

# Validate heating rates
python tests/radiation_profiles.py --case tropical_atmosphere
```

### Benchmark Cases
- **Clear-sky**: Known analytical solutions
- **Radiative-convective**: Equilibrium temperature profiles
- **Cloud effects**: Impact of cloud optical properties

### Expected Results

#### Clear-Sky Cooling
- **Tropics**: Net cooling ~1-2 K/day
- **Mid-latitudes**: Seasonal variation in cooling rates
- **Polar regions**: Strong LW cooling, weak SW heating

#### Diurnal Cycle
- **Surface temperature**: ~10-20K diurnal range
- **Boundary layer**: Deepening during day, shallowing at night
- **Cloud formation**: Radiative cooling drives nocturnal clouds

## Performance Characteristics

### Computational Cost
- **Simple grey**: O(N_grid) per radiation call
- **RRTMG**: O(N_grid × N_bands × N_streams)
- **Frequency**: Every 5-15 minutes of simulation time

### Memory Requirements
- **Flux arrays**: 3D arrays for up/downward fluxes
- **Optical properties**: 3D arrays for absorption coefficients
- **Spectral data**: Band-specific arrays (RRTMG)

## Scientific References

### Radiative Transfer Theory
- Chandrasekhar (1960) - Radiative Transfer
- Goody & Yung (1989) - Atmospheric Radiation
- Liou (2002) - An Introduction to Atmospheric Radiation

### Parameterization Methods
- Fels & Schwarzkopf (1975) - Grey atmosphere approximation
- Ritter & Geleyn (1992) - Simple radiation schemes
- Mlawer et al. (1997) - RRTM development

### GCM Applications
- Kiehl & Trenberth (1997) - Earth's radiation budget
- Collins et al. (2006) - CAM radiation schemes
- Iacono et al. (2008) - RRTMG validation

This module provides the radiative forcing needed for realistic atmospheric temperature profiles and diurnal cycles.
