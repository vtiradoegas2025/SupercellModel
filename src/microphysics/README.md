# Cloud Microphysics Module

This module implements parameterizations for cloud and precipitation microphysics, handling phase changes, particle growth, and precipitation formation in atmospheric simulations.

## Overview

Microphysics schemes simulate the evolution of hydrometeor species through:
- **Nucleation**: Formation of cloud droplets/ice crystals
- **Growth processes**: Condensation, deposition, riming, aggregation
- **Phase changes**: Freezing, melting, sublimation
- **Precipitation**: Fallout and sedimentation

## Architecture

```
src/microphysics/
├── microphysics.cpp           # Module coordinator
├── factory.cpp/.hpp           # Scheme factory
├── base/
│   └── thermodynamics.cpp/.hpp # Phase change thermodynamics
└── schemes/                   # Microphysics implementations
    ├── kessler/              # Warm-rain scheme
    ├── lin/                  # Ice microphysics
    ├── thompson/             # Aerosol-aware mixed-phase
    └── milbrandt/            # Double-moment scheme
```

## Hydrometeor Species

### Prognostic Variables
- **q_v**: Water vapor mixing ratio (kg/kg)
- **q_c**: Cloud water mixing ratio (kg/kg)
- **q_r**: Rain water mixing ratio (kg/kg)
- **q_i**: Cloud ice mixing ratio (kg/kg)
- **q_s**: Snow mixing ratio (kg/kg)
- **q_g**: Graupel mixing ratio (kg/kg)
- **q_h**: Hail mixing ratio (kg/kg)

### Optional Number Concentrations (Double-Moment)
- **N_c**: Cloud droplet number concentration (m⁻³)
- **N_r**: Rain drop number concentration (m⁻³)
- **N_i**: Ice crystal number concentration (m⁻³)
- **N_s**: Snow aggregate number concentration (m⁻³)
- **N_g**: Graupel number concentration (m⁻³)
- **N_h**: Hailstone number concentration (m⁻³)

## Scheme Implementations

### Kessler Warm-Rain Scheme

**Species**: q_v, q_c, q_r (3-moment: adds N_r)

**Key Processes:**
```cpp
// Autoconversion: cloud → rain
dq_r/dt = c_auto * q_c * (q_c / q_c0)^α when q_c > q_c0

// Accretion: rain collects cloud
dq_r/dt = c_accr * q_r * q_c

// Evaporation: rain evaporates in subsaturated air
dq_r/dt = c_evap * q_r * (q_vs - q_v) / q_vs
```

**Parameters:**
- q_c0 = 1.0×10⁻³ kg/kg (autoconversion threshold)
- c_auto = 1.0×10⁻³ s⁻¹ (autoconversion rate)
- c_accr = 2.2 s⁻¹ (accretion rate)
- c_evap = 5.0×10⁻³ s⁻¹ (evaporation rate)

### Lin et al. Ice Microphysics

**Species**: q_v, q_c, q_r, q_i, q_s (7 variables)

**Ice Processes:**
- **Homogeneous freezing**: Water → ice at T < T_hom
- **Heterogeneous nucleation**: Ice formation on aerosols
- **Deposition/sublimation**: Vapor ↔ ice phase changes
- **Aggregation**: Ice particle collisions → snow
- **Riming**: Ice particles collect supercooled water
- **Melting**: Ice → water at T > 0°C

**Size Distributions:**
- Cloud droplets: Gamma distribution
- Rain: Marshall-Palmer (exponential)
- Ice: Exponential with temperature-dependent slopes

### Thompson Aerosol-Aware Scheme

**Species**: q_v, q_c, q_r, q_i, q_s, q_g (6 variables)

**Aerosol Coupling:**
- **CCN activation**: Cloud droplet formation from aerosols
- **Ice nucleation**: Deposition/immersion freezing on aerosols
- **Size-dependent growth**: Aerosol influence on particle sizes

**Key Features:**
- **Aerosol-aware CCN**: Meyers et al. (1992) parameterization
- **Temperature-dependent nucleation**: Cooper (1986) for ice
- **Variable density**: Different ice habits (dendritic, plates, columns)

**Size Distributions:**
- Cloud droplets: Gamma (shape parameter from aerosols)
- Precipitation: Gamma distributions with prognostic intercepts

### Milbrandt-Yau Double-Moment Scheme

**Species**: q_v, q_c, q_r, q_i, q_s, q_g, q_h (7 variables)
**Number concentrations**: N_c, N_r, N_i, N_s, N_g, N_h (6 additional)

**Advanced Features:**
- **Prognostic number concentrations**: Full DSD evolution
- **Shape parameters**: Gamma distribution shape evolution
- **Multi-category ice**: Separate treatment of ice habits
- **Hail prediction**: Separate hail category with size sorting

## Thermodynamics Implementation

### Phase Change Calculations

#### Saturation Vapor Pressures
```cpp
// Tetens formula for water vapor
e_s(T) = e_s0 * exp(a * (T - T0) / (T - b))

// Ice saturation (below 0°C)
e_si(T) = e_s0 * exp(a_i * (T - T0) / (T - b_i))
```

#### Latent Heating
```cpp
// Condensation/evaporation
dθ/dt = (L_v / c_p) * (dq_v/dt) * (1 + (R_v q_v)/(R_d (1 - q_t)))

// Freezing/melting
dθ/dt = (L_f / c_p) * (dq_ice/dt)

// Deposition/sublimation
dθ/dt = (L_s / c_p) * (dq_ice/dt)
```

Where:
- L_v = 2.5×10⁶ J/kg (latent heat of vaporization)
- L_f = 3.34×10⁵ J/kg (latent heat of fusion)
- L_s = L_v + L_f (latent heat of sublimation)

## Numerical Implementation

### Time Integration
Microphysics tendencies are computed diagnostically and applied with the dynamics time step:

```cpp
// In dynamics loop
step_dynamics(dt);  // Apply dynamics tendencies
step_microphysics(dt);  // Compute and apply microphysics
```

### Stability Considerations
- **Positive definiteness**: Ensure q ≥ 0 after each process
- **Conservation**: Water mass conservation across phase changes
- **Numerical diffusion**: Implicit treatment of fast processes

### Process Coupling
Processes are applied in order of increasing time scales:
1. **Fast**: Condensation/deposition (τ ~ seconds)
2. **Medium**: Nucleation, freezing (τ ~ 10²-10³ s)
3. **Slow**: Collection, fallout (τ ~ 10³-10⁴ s)

## Configuration

### Scheme Selection
```yaml
microphysics:
  scheme: "thompson"     # kessler, lin, thompson, milbrandt
  dt_micro: 30.0         # Microphysics time step (s)
```

### Scheme-Specific Parameters
```yaml
# Kessler
microphysics:
  kessler:
    qc0: 0.001           # Autoconversion threshold
    c_auto: 0.001        # Autoconversion rate
    c_accr: 2.2          # Accretion rate
    c_evap: 0.003        # Evaporation rate

# Thompson
microphysics:
  thompson:
    nccn: 100.0e6        # CCN concentration (m⁻³)
    nino: 0.1e6          # Ice nuclei concentration (m⁻³)
```

## Validation and Testing

### Benchmark Tests
```bash
# Test microphysics schemes
python tests/validate_microphysics.py --scheme thompson --case warm_rain

# Compare schemes on same initial condition
python tests/compare_schemes.py --schemes kessler thompson --case supercell
```

### Physical Consistency Checks
- **Water conservation**: ∫(q_v + q_c + q_r + q_i + q_s + q_g + q_h) conserved
- **Energy conservation**: Latent heating balances phase changes
- **Precipitation rates**: Consistent with observed climatology
- **Size distributions**: Realistic DSD parameters

## Performance Characteristics

### Computational Cost
- **Kessler**: ~5-10% of total simulation time
- **Lin**: ~10-15% of total simulation time
- **Thompson**: ~15-20% of total simulation time
- **Milbrandt**: ~25-30% of total simulation time

### Memory Requirements
- **Single-moment**: 5-7 3D arrays
- **Double-moment**: 11-13 3D arrays
- **Cache efficiency**: Good locality for prognostic variables

## Scientific References

### Kessler (1969)
- Original warm-rain parameterization
- Simple, computationally efficient
- Widely used in idealized studies

### Lin et al. (1983)
- First comprehensive bulk ice scheme
- Temperature-dependent processes
- Foundation for modern microphysics

### Thompson et al. (2008)
- Aerosol-aware microphysics
- NCAR operational use
- Advanced ice particle properties

### Milbrandt & Yau (2005)
- Double-moment approach
- Prognostic number concentrations
- High-fidelity precipitation simulation

This module provides the foundation for cloud and precipitation processes in atmospheric simulations.
