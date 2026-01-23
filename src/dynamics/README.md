# Atmospheric Dynamics Module

This module implements the compressible non-hydrostatic Euler equations for atmospheric dynamics, providing the core fluid motion solver for convective storm simulations.

## Overview

The dynamics module solves the fundamental equations governing atmospheric motion:
- **Mass continuity**: Conservation of air mass
- **Momentum equations**: Newton's laws for fluid motion
- **Thermodynamic equation**: Conservation of potential temperature
- **Equation of state**: Ideal gas relationship

## Governing Equations

### Compressible Euler Equations (Non-Hydrostatic)

In cylindrical coordinates (r, θ, z), the equations are:

**Continuity (Mass conservation):**
```
∂ρ/∂t + (1/r) ∂(r ρ u)/∂r + (1/r) ∂(ρ v)/∂θ + ∂(ρ w)/∂z = 0
```

**Radial momentum:**
```
∂u/∂t + u ∂u/∂r + (v/r) ∂u/∂θ + w ∂u/∂z - v²/r = - (1/ρ) ∂p/∂r + ν ∇²u + F_u
```

**Azimuthal momentum:**
```
∂v/∂t + u ∂v/∂r + (v/r) ∂v/∂θ + w ∂v/∂z + (u v)/r = - (1/(ρ r)) ∂p/∂θ + ν ∇²v + F_v
```

**Vertical momentum:**
```
∂w/∂t + u ∂w/∂r + (v/r) ∂w/∂θ + w ∂w/∂z = - (1/ρ) ∂p/∂z - g + ν ∇²w + F_w + B
```

**Potential temperature:**
```
∂θ/∂t + u ∂θ/∂r + (v/r) ∂θ/∂θ + w ∂θ/∂z = ν ∇²θ + Q_radiation + Q_microphysics
```

**Equation of state:**
```
p = ρ R_d T_v = ρ R_d θ (π / p_00)^{R_d/c_p}
```

Where:
- ρ: density, u,v,w: velocity components, p: pressure, θ: potential temperature
- g: gravity, ν: kinematic viscosity, B: buoyancy, Q: source terms
- R_d: dry air gas constant, c_p: specific heat at constant pressure

## Architecture

```
src/dynamics/
├── dynamics.cpp              # Module coordinator
├── factory.cpp/.hpp          # Dynamics scheme factory
└── schemes/                  # Dynamics implementations
    ├── supercell/            # Supercell-specific dynamics
    └── tornado/              # Tornado-scale dynamics
```

## Scheme Implementations

### Supercell Dynamics Configuration

**Target applications:** Classic supercell thunderstorms, mesoscale convective systems

**Features:**
- Standard compressible formulation
- Optimized for convection-permitting resolutions (Δx ≈ 1km)
- Balanced for mid-latitude severe weather

**Configuration:**
```yaml
dynamics:
  scheme: "supercell"
  compressibility: true      # Use compressible equations
  hydrostatic_balance: false # Non-hydrostatic formulation
  coriolis: false           # Neglect Coriolis for idealized cases
```

### Tornado Dynamics Configuration

**Target applications:** Tornado-scale vortices, high-resolution simulations

**Features:**
- Enhanced numerical stability for fine grids
- Specialized for Δx < 100m resolutions
- Optimized for strong vertical gradients

**Configuration:**
```yaml
dynamics:
  scheme: "tornado"
  compressibility: true
  high_order_advection: true  # Higher-order spatial accuracy
  implicit_diffusion: true    # Better stability for small Δx
```

## Numerical Implementation

### Spatial Discretization

**Finite-difference operators:**
```cpp
// Momentum advection (2nd-order centered)
A_u = u ∂u/∂r + (v/r) ∂u/∂θ + w ∂u/∂z - v²/r

// Pressure gradient
G_r = - (1/ρ) ∂p/∂r

// Buoyancy term
B = g (θ_v' / θ_{v0})  // Virtual potential temperature perturbation
```

**Conservative form:**
```cpp
∂(ρ u)/∂t + ∇ · (ρ u ⊗ u) + ∇p + ∇ · τ = ρ g k̂ + ρ B
```

### Time Integration

**Split-explicit approach:**
1. **Acoustic modes**: Implicit or semi-implicit treatment
2. **Advection modes**: Explicit high-order schemes
3. **Gravity modes**: Explicit with CFL-limited time steps

**RK3 integration for advection:**
```cpp
u^{(1)} = u^n + Δt L_{adv}(u^n)
u^{(2)} = (3/4)u^n + (1/4)u^{(1)} + (1/4)Δt L_{adv}(u^{(1)})
u^{n+1} = (1/3)u^n + (2/3)u^{(2)} + (2/3)Δt L_{adv}(u^{(2)})
```

### Pressure Solver

**Poisson equation for pressure:**
```cpp
∇²π' = (g/θ_0) ∂θ_v/∂z + ∇ · (∇ · (u ⊗ u)) / c_s²
```

Where π' is Exner pressure perturbation, c_s is sound speed.

**Solution methods:**
- **Direct solver**: For small domains
- **Iterative solver**: Conjugate gradient for large domains
- **Multigrid**: For optimal scaling

### Boundary Conditions

**Lateral boundaries:**
```cpp
// Radiative condition for waves
∂ψ/∂t + c ∂ψ/∂n = 0  // c = phase speed

// Open boundary for scalars
∂θ/∂n = 0, ∂q/∂n = 0

// Reflective for normal velocity
u_n = 0
```

**Top boundary (Rayleigh damping):**
```cpp
// Sponge layer relaxation
∂ψ/∂t = - (ψ - ψ_ref) / τ(z)
τ(z) = τ_0 exp(-((z - z_damp)/H_damp)²)
```

**Bottom boundary:**
```cpp
// Free-slip condition
w = 0
∂u/∂z = 0, ∂v/∂z = 0

// Surface fluxes from PBL scheme
τ_{xz} = -ρ u_*², τ_{yz} = -ρ v_*²
```

## Physical Parameterizations

### Buoyancy

**Virtual potential temperature:**
```cpp
θ_v = θ (1 + ε q_v - q_c - q_r - q_i - q_s - q_g - q_h)
```

Where ε = (R_v / R_d) - 1 ≈ 0.608

**Buoyancy acceleration:**
```cpp
B = g (θ_v' / θ_{v0} - q_l)  // q_l = liquid/ice condensate
```

### Coriolis Force (Optional)

**For large domains:**
```cpp
F_u = f v, F_v = -f u  // f = 2Ω sinφ
```

### Acoustic Filtering

**Time-splitting for stability:**
- **Slow modes**: Explicit advection of momentum/thermo
- **Fast modes**: Implicit acoustic integration

## Configuration Parameters

### Basic Settings
```yaml
dynamics:
  scheme: "supercell"
  dt: 0.1                # Dynamics time step (s)
  cfl_max: 0.8           # Maximum CFL number
  compressibility: true  # Use compressible equations
  coriolis: false        # Include Coriolis force
```

### Advanced Settings
```yaml
dynamics:
  # Pressure solver
  pressure_solver: "iterative"  # direct, iterative, multigrid
  max_iterations: 100
  tolerance: 1e-6

  # Boundary conditions
  lateral_bc: "radiative"       # radiative, periodic, outflow
  top_damping: true
  damping_height: 10000.0       # Damping layer height (m)
  damping_strength: 300.0       # Relaxation time (s)

  # Numerical options
  advection_order: 5            # 3 or 5 for WENO
  implicit_acoustics: true      # Implicit acoustic solver
```

## Implementation Details

### Memory Management

**3D field allocation:**
```cpp
std::vector<std::vector<std::vector<float>>> u(NR,
    std::vector<std::vector<float>>(NTH,
        std::vector<float>(NZ, 0.0f)));
```

**Cache optimization:**
- Contiguous memory access patterns
- SOA (struct-of-arrays) for vectorization
- Minimal memory allocations during runtime

### Parallel Scaling

**Domain decomposition:**
- Natural splitting in radial direction
- Periodic boundaries in azimuthal direction
- Surface/top boundary handling in vertical

**Communication patterns:**
- Ghost cell exchanges for finite differences
- All-reduce for global diagnostics
- Load balancing for varying resolution

### Performance Optimization

**Computational kernels:**
- Vectorized inner loops with SIMD
- Cache-blocking for large arrays
- Precomputed geometric factors

**Scalability:**
- O(NR × NTH × NZ) scaling for 3D operations
- O(NR × NTH × NZ × log NZ) for pressure solver
- Memory bandwidth limited for large domains

## Validation and Testing

### Benchmark Cases
```bash
# Test dynamics schemes
python tests/test_dynamics.py --scheme supercell --resolution 256x128x128

# Validate conservation properties
python tests/conservation_test.py --case thermal_bubble
```

### Convergence Studies
- **Grid refinement**: Error reduction with Δx
- **Time step convergence**: RK3 accuracy verification
- **Conservation**: Mass, momentum, energy invariants

### Expected Performance
- **Small domains (32³)**: >1000 time steps/second
- **Large domains (256×128×128)**: ~10-50 time steps/second
- **Memory usage**: ~4GB for 256×128×128 single precision

## Scientific References

### Fundamental Equations
- Klemp & Wilhelmson (1978) - Compressible storm dynamics
- Rotunno & Klemp (1985) - Supercell rotation and propagation
- Wicker & Skamarock (2002) - Time-splitting for elastic models

### Numerical Methods
- Durran (1999) - Numerical Methods for Wave Equations
- Skamarock & Klemp (2008) - Terrain-following coordinates
- Weller et al. (2013) - Conservative finite-difference schemes

### Validation Studies
- Bryan & Fritsch (2002) - Benchmark moist nonhydrostatic models
- Bryan et al. (2003) - Resolution requirements for deep convection
- Klemp et al. (2007) - Upper gravity wave absorbing layer

This module provides the fundamental fluid dynamics solver for atmospheric simulations.
