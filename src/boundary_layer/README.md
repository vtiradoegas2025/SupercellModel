# Planetary Boundary Layer Module

This module implements parameterizations of the planetary boundary layer (PBL), which couples the atmosphere to the Earth's surface through turbulent fluxes of momentum, heat, and moisture.

## Overview

The PBL is the lowest 1-3 km of the atmosphere where:
- **Surface drag**: Momentum transfer from wind to surface
- **Heat exchange**: Sensible and latent heat fluxes
- **Moisture transfer**: Evaporation and condensation
- **Turbulent mixing**: Vertical transport of scalars and momentum

## Architecture

```
src/boundary_layer/
├── boundary_layer.cpp       # Module coordinator
├── factory.cpp/.hpp         # PBL scheme factory
├── base/
│   └── surface_fluxes.cpp/.hpp # Surface flux calculations
└── schemes/                 # PBL implementations
    ├── slab/                # Simple mixed-layer model
    ├── ysu/                 # Yonsei University scheme
    └── mynn/                # Mellor-Yamada-Nakanishi-Niino
```

## PBL Physics Fundamentals

### Turbulent Fluxes

**Reynolds stresses:**
```
-ρ <u' w'> = ρ u_*²  (momentum flux)
```

**Scalar fluxes:**
```
-ρ c_p <θ' w'> = ρ c_p θ_*²  (heat flux)
-ρ L_v <q' w'> = ρ L_v q_*²  (moisture flux)
```

Where u_*, θ_*, q_* are surface scales.

### Boundary Layer Height

**Mixed layer depth:**
```cpp
h = √(2 ∫_0^∞ (θ_v - θ_{v,sfc}) / (dθ_v/dz|_top) dz)
```

**Entrainment rate:**
```cpp
dh/dt = w_e = (w_*³ / (5 h Δb))^{1/3}
```

Where w_* is convective velocity scale.

## Scheme Implementations

### Slab Model

**Simple mixed-layer approach:**
```cpp
// Constant flux layer
∂u/∂t = (u_*² / h) × stability_function
∂θ/∂t = (θ_*² / h) × stability_function
∂q/∂t = (q_*² / h) × stability_function
```

**Configuration:**
```yaml
boundary_layer:
  scheme: "slab"
  pbl_height_m: 1000.0       # Fixed boundary layer height
```

**Features:**
- **Simple**: Minimal computational cost
- **Stable**: No convergence issues
- **Limited**: No vertical structure

### Yonsei University (YSU) Scheme

**Non-local mixing with explicit entrainment:**

**Key features:**
- **Counter-gradient fluxes**: Organized transport in convective conditions
- **Explicit entrainment**: Cloud-top radiative cooling effects
- **Large-eddy simulation**: Inspired mixing length formulation

**Turbulent diffusivity:**
```cpp
K_m = κ u_* z × Φ_m(ζ) × (1 - z/h)^2
```

Where Φ_m is the stability function.

**Configuration:**
```yaml
boundary_layer:
  scheme: "ysu"
  entrainment: true          # Include entrainment effects
  cloud_top_cooling: true    # Cloud-top radiative effects
```

### Mellor-Yamada-Nakanishi-Niino (MYNN)

**TKE-based local closure:**

**Prognostic TKE equation:**
```cpp
d(e)/dt = P + B - ε + transport
```

Where:
- **P**: Shear production (u_i' u_j' ∂U_i/∂x_j)
- **B**: Buoyancy production (β θ' w')
- **ε**: Dissipation (c_ε e^{3/2}/l)

**Eddy diffusivity:**
```cpp
K_m = c_μ √e l_m
K_h = K_m / Pr_t
```

**Configuration:**
```yaml
boundary_layer:
  scheme: "mynn"
  tke_min: 1e-6             # Minimum TKE value
  length_scale: "mynn"      # MYNN length scale formulation
  cloud_pdf: "single"       # Cloud PDF for moist processes
```

## Surface Layer Parameterization

### Monin-Obukhov Similarity Theory (MOST)

**Universal functions for stable/unstable regimes:**

**Momentum:**
```cpp
φ_m(ζ) = 1 + β ζ    (stable, ζ > 0)
φ_m(ζ) = (1 - γ ζ)^{-1/4}  (unstable, ζ < 0)
```

**Heat:**
```cpp
φ_h(ζ) = Pr_t φ_m(ζ)  (approximately)
```

Where ζ = z/L, L is Obukhov length.

### Surface Flux Calculation

**Iterative solution for surface scales:**
```cpp
// Initial guess
u_* = κ U₁ / ln(z₁/z₀)
θ_* = κ (θ₁ - θ_sfc) / ln(z₁/z₀)
L = - (θ_* u_*^2) / (κ g θ_v)

// Iterate with stability corrections
u_* = κ U₁ / (ln(z₁/z₀) - ψ_m(ζ))
θ_* = κ (θ₁ - θ_sfc) / (ln(z₁/z₀) - ψ_h(ζ))
```

## Implementation Details

### Vertical Diffusion

**Implicit solution for stability:**
```cpp
# Crank-Nicolson scheme
(A - (Δt/2) K ∇²) φ^{n+1} = (A + (Δt/2) K ∇²) φ^n + S
```

Where A is advection operator, K is diffusivity, S is source terms.

### Non-Local Transport (YSU)

**Counter-gradient term:**
```cpp
∂θ/∂z|_cgt = γ_cgt (w_θ_surf / w_*)
```

Where γ_cgt is the counter-gradient coefficient.

### Cloud-Top Entrainment

**Radiative cooling effect:**
```cpp
dθ/dt|_cloudtop = - (F_rad / (ρ c_p)) * (∂θ/∂z)^{-1}
```

## Configuration

### Basic PBL Setup
```yaml
boundary_layer:
  scheme: "ysu"             # slab, ysu, mynn
  dt_pbl: 60.0              # PBL time step (s)
  apply_surface_fluxes: true
```

### Surface Parameters
```yaml
surface:
  z0m: 0.01                 # Momentum roughness length (m)
  z0h: 0.001                # Heat roughness length (m)
  albedo: 0.2               # Surface albedo
  emissivity: 0.98          # Surface emissivity
```

### Advanced Options
```yaml
boundary_layer:
  # YSU specific
  entrainment_coefficient: 0.2
  cloud_top_radiation: true

  # MYNN specific
  tke_advection: true       # Include TKE advection
  length_scale_damping: true # Damp length scale near surface
```

## Validation and Testing

### Surface Layer Tests
```bash
# Test MOST functions
python tests/test_most.py --stability neutral

# Validate surface fluxes
python tests/test_fluxes.py --wind_speed 10 --temperature_diff 2
```

### PBL Profile Tests
```bash
# Test mixed layer development
python tests/test_pbl.py --scheme ysu --forcing convective

# Compare schemes
python tests/compare_pbl.py --schemes slab ysu mynn --case neutral
```

### Expected Results

#### Neutral Boundary Layer
- **Wind profile**: Logarithmic from surface to boundary layer top
- **Temperature**: Constant potential temperature
- **Fluxes**: Constant with height

#### Convective Boundary Layer
- **Deep mixing**: Well-mixed profiles up to entrainment zone
- **Cloud formation**: Cumulus clouds with cold pools
- **Entrainment**: Cloud-top mixing with free troposphere

#### Stable Boundary Layer
- **Weak mixing**: Sharp gradients near surface
- **Decoupling**: Surface inversion layer
- **Gravity waves**: Possible wave generation

## Performance Characteristics

### Computational Cost
- **Slab**: O(N_grid) per PBL step
- **YSU**: O(N_grid × N_levels) with vertical iteration
- **MYNN**: O(N_grid × N_levels) with TKE solver

### Memory Requirements
- **Slab**: Minimal additional storage
- **YSU**: Working arrays for vertical integration
- **MYNN**: Additional TKE prognostic variable

## Scientific References

### PBL Theory
- Stull (1988) - An Introduction to Boundary Layer Meteorology
- Garratt (1992) - The Atmospheric Boundary Layer
- Wyngaard (2010) - Turbulence in the Atmosphere

### Parameterization Development
- Blackadar (1962) - The vertical distribution of wind and turbulent exchange
- Mellor & Yamada (1982) - Development of a turbulence closure model
- Hong et al. (2006) - A new vertical diffusion package

### MOST and Surface Layer
- Monin & Obukhov (1954) - Basic laws of turbulent mixing
- Businger et al. (1971) - Flux-profile relationships in the atmospheric surface layer
- Dyer (1974) - A review of flux-profile relationships

This module provides the critical coupling between atmospheric dynamics and surface processes.
