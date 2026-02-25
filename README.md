# SupercellModel

[![C++](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](#license)

A high-performance, research-grade atmospheric simulation framework for **supercell and tornadic storm research**. SupercellModel implements compressible, non-hydrostatic equations in cylindrical coordinates with a modular physics architecture (**C++17**) and an active **native Vulkan rendering path** for local analysis.

> **Status:** Work in progress / research prototype. APIs, file formats, and results may change.

See **[docs/STATUS.md](docs/STATUS.md)**
---

## Project Goal (CM1-lite on Personal Hardware)

SupercellModel is an attempt to build a **CM1-lite style research model** that can run on workstation/personal hardware while keeping the physics architecture modular and testable.

This project is not trying to replace CM1; CM1 remains the reference-standard system. The goal here is a pragmatic subset with strong developer ergonomics and reproducible validation workflows.

### Why this codebase is C++/modern-tooling oriented

Fortran remains strong for production atmospheric modeling, but this project emphasizes:
- explicit memory-layout control for large 3D fields (`Field3D` flattened storage)
- easier integration with modern profiling/testing/tooling workflows
- modular factory-based scheme experimentation across physics components
- direct native rendering path integration (`vulkan/`) for local analysis loops

### Scope for this repo

- Runtime simulation path: `bin/tornado_sim`
- Native viewer path: `bin/vulkan_viewer`
- Strict exported-field validation contract and guard tooling
- Modular physics scheme comparison for storm-scale research

### Non-goals (current)

- Full CM1 diagnostic breadth parity in one release
- Operational forecasting-grade completeness
- Eliminating all scientific calibration work (terrain/chaos tuning is ongoing)

### Current Coverage Snapshot (February 25, 2026)

- CM1-style field contract: `99` fields tracked
- Exported now: `50`
- Backlog/not implemented: `49`
- Required-now coverage: `20/20` exported
- Validation gates currently used:
  - `make test`
  - `make test-backend-physics`
  - `make test-soundings`

For the current gap audit and roadmap, see [docs/STATUS.md](docs/STATUS.md).

---

## Key Features

#### Core Simulation Engine (C++17)
- **Compressible Dynamics**: Full Euler equations without hydrostatic approximation
- **Cylindrical Grid**: Optimized for storm-scale simulations (r, θ, z coordinates)
- **RK3 Time Integration**: 3rd-order Runge-Kutta with CFL-limited adaptive stepping
- **Modular Physics**: Factory-based architecture for easy parameterization comparison
- **Performance Optimizations**: 
  - Flattened array storage (Field3D) for improved cache locality
  - OpenMP multi-core parallelization
  - SIMD-ready architecture for vectorized operations
  - CPU-specific compiler optimizations

#### Physics Parameterizations
- **Microphysics**: 4 schemes (Kessler warm-rain, Thompson, Lin, Milbrandt)
- **Boundary Layer**: 3 schemes (YSU, MYNN, slab)
- **Turbulence**: 2 schemes (Smagorinsky-Lilly, TKE prognostic)
- **Radiation**: 1 scheme implemented (`simple_grey`; RRTMG planned)
- **Radar Operators**: Forward simulators for Z/Ze, V_r, Z_DR, K_DP, ρ_HV

#### 3D Visualization (Native Vulkan)
- **Vulkan Viewer (`vulkan/`)**: Native Vulkan renderer for direct NPY ingest from backend exports
- **Volume Rendering**: GPU-accelerated ray marching through atmospheric volumes
- **Coordinate Transformation**: Cylindrical simulation data → Cartesian visualization
- **Interactive Viewer**: Real-time exploration via native orbital camera controls

### Example Output
Simulated classic supercell thunderstorm showing:
- Rotating updraft with vertical vorticity > 0.01 s⁻¹
- Asymmetric precipitation with hook echo formation
- Radar reflectivity > 50 dBZ in convective core
- Doppler velocity dipole structure

---

## System Architecture

### Governing Equations

The model solves the compressible Euler equations in cylindrical coordinates:

```
∂u/∂t + ∇·(u⊗u) + (1/ρ₀)∇p' + gθ'/θ₀ k̂ = -∇·τ + F_buoyancy
∂w/∂t + ∇·(u⊗w) + (1/ρ₀)∇p' + g(θ'/θ₀ - q_t) k̂ = -∇·τ
∂θ/∂t + ∇·(uθ) = Q_radiation + Q_microphysics + ∇·(K_θ ∇θ)
∂q_v/∂t + ∇·(uq_v) = -C - E + ∇·(K_q ∇q_v)
```

Where ρ₀ is base-state density, p' is perturbation pressure, τ is sub-grid stress tensor, and Q terms represent source/sink processes.

### Coordinate System & Grid
- **Primary Grid**: Cylindrical (r, θ, z) with uniform spacing
- **Visualization**: Cartesian transformation for 3D rendering
- **Resolution**: Configurable NR × NTH × NZ (typically 128×128×64 to 256×256×128)
- **Domain**: r ∈ [0, r_max], θ ∈ [0, 2π], z ∈ [0, z_max]

### Boundary Conditions
- **Lateral (r = r_max)**: Radiative outflow with open boundary
- **Top (z = z_max)**: Rayleigh damping sponge layer
- **Bottom (z = 0)**: Free-slip with Monin-Obukhov surface fluxes

---

## Repository Structure

```
SupercellModel/
├── src/                    # Core simulation code (C++17)
│   ├── core/              # Runtime coordinators and main executable
│   ├── dynamics/          # Dynamics schemes
│   ├── microphysics/      # Cloud physics parameterizations
│   ├── boundary_layer/    # Planetary boundary layer schemes
│   ├── turbulence/        # Sub-grid turbulence closures
│   ├── radiation/         # Radiative transfer schemes
│   ├── radar/            # Forward radar operators
│   └── numerics/         # Advection/diffusion/time-stepping
├── include/               # Public headers and API
├── configs/               # Pre-configured simulation setups
├── vulkan/                # Native Vulkan renderer
├── docs/                  # Technical documentation
├── bin/                   # Build outputs
└── data/                  # Sample data and exports
```

---

## Prerequisites & Installation

### System Requirements
- **C++ Compiler**: C++17 compliant (clang++ ≥ 9.0, g++ ≥ 7.0)
- **OpenMP**: Optional but recommended for multi-core performance (libomp on macOS via Homebrew)
- **Graphics**: Vulkan-capable GPU/driver stack (MoltenVK on macOS)
- **Memory**: 8GB+ RAM recommended for production simulations
- **CPU**: Multi-core processor recommended (4+ cores for optimal performance)

### Dependencies
- **C++**: Standard Template Library, OpenMP (optional)
- **Vulkan path**: vulkan-loader, Vulkan headers, MoltenVK, GLFW, glslang

### Build & Install

Build/run commands are consolidated in **Quick Start & Command Reference** near the bottom of this README.
OpenMP is auto-detected; on macOS you can install `libomp` via `brew install libomp`.

---

## Example Configuration Snippets

**Classic Supercell** (Weisman-Klemp):
```yaml
microphysics:
  scheme: thompson
boundary_layer:
  scheme: ysu
turbulence:
  scheme: smagorinsky
  dt_sgs: 1.0
  Cs: 0.18
  Pr_t: 0.7
  Sc_t: 0.7
environment:
  cape_target_jkg: 2500
```

**High-Resolution Tornado Case**:
```yaml
grid:
  nr: 512
  nth: 256
  nz: 256
microphysics:
  scheme: milbrandt
radiation:
  scheme: simple_grey
```

---

## Physics Modules

### Microphysics Parameterizations
| Scheme | Species | Features | Use Case |
|--------|---------|----------|----------|
| **Kessler** | q_v, q_c, q_r (+ q_g, q_h in this implementation) | Warm-rain core, extended mixed-phase pathways | Idealized studies |
| **Thompson** | q_v, q_c, q_r, q_i, q_s, q_g, q_h (+ internal N_i) | Mixed-phase with internal ice-number evolution | Realistic storms |
| **Lin** | q_v, q_c, q_r, q_i, q_s, q_g, q_h | Bulk mixed-phase with graupel/hail categories | Mixed-phase research |
| **Milbrandt-Yau** | q_v, q_c, q_r, q_i, q_s, q_g, q_h + internal N_r/N_i/N_s/N_g/N_h | Double-moment style | Size distribution studies |

### Boundary Layer Schemes
- **YSU**: Non-local mixing, explicit entrainment, counter-gradient fluxes
- **MYNN**: TKE-based local closure, stability-dependent length scales
- **Slab**: Simple mixed-layer model for idealized simulations

### Turbulence Closures
- **Smagorinsky-Lilly**: Eddy-viscosity model (C_s = 0.18)
- **TKE Prognostic**: Full prognostic equation with production/dissipation

### Radar Forward Operators
- **Reflectivity (Z/Ze)**: Rayleigh scattering from hydrometeors
- **Doppler Velocity (V_r)**: Radial component of wind + fall speeds
- **Polarimetric (Z_DR, K_DP, ρ_HV)**: Dual-polarization variables

---

## Validation & Performance

### Physical Validation
- **Mass Conservation**: <0.1% drift over 2-hour simulations
- **Energy Conservation**: <5% total energy change in stable cases
- **Storm Structure**: Realistic supercell morphology vs. literature
- **Radar Signatures**: Z-V_r relationships match theoretical expectations

### Regression Tests
Primary command workflows are consolidated in **Quick Start & Command Reference** below.

### Performance Benchmarks
- **Small Test Grid** (64×64×32): ~10,000-15,000 time steps/hour (with optimizations)
- **Production Grid** (256×128×128): ~100-150 time steps/hour (with optimizations)
- **Memory Usage**: ~8GB for large domains (optimized with flattened array storage)
- **Visualization**: 30-60 FPS interactive rendering
- **Optimizations**: 
  - Compiler: `-O3 -march=native -mtune=native` for CPU-specific optimizations
  - OpenMP: Multi-core parallelization (4-8x speedup on multi-core systems)
  - Memory: Flattened array storage reduces allocation overhead by ~99.9%
  - Cache: Contiguous memory layout improves cache locality (1.5-3x speedup)

### Known Limitations
- Terrain module is integrated, but broader physics validation/calibration is ongoing
- Chaos/stochastic perturbations are integrated; ensemble calibration and tuning are ongoing
- Some advanced radar operators in development
- APIs and file formats may change

---

## Documentation & Resources

### Documentation Structure
- **[Technical Reference](docs/README.md)** - Complete system documentation
- **[Scientific Foundation](docs/foundationalScience.md)** - Literature references
- **[API Reference](include/README.md)** - Code interfaces and headers
- **[Source Code Guide](src/README.md)** - Implementation details

### Module-Specific Documentation
- **[Microphysics](src/microphysics/README.md)** - Cloud physics schemes
- **[Boundary Layer](src/boundary_layer/README.md)** - PBL parameterizations
- **[Turbulence](src/turbulence/README.md)** - SGS closures
- **[Radar](src/radar/README.md)** - Forward operators
- **[Numerics](src/numerics/README.md)** - Numerical methods
- **[Radiation](src/radiation/README.md)** - Radiative transfer
- **[Chaos](src/chaos/README.md)** - Stochastic perturbations
- **[Terrain](src/terrain/README.md)** - Orographic effects

### Development & Testing
- **[Project Status](docs/STATUS.md)** - Active gap audit and current validation state
- **[Vulkan Viewer Guide](vulkan/README.md)** - Native rendering path and runtime options

---

## Quick Start & Command Reference

### Clone and Build

```bash
git clone https://github.com/vtiradoegas2025/SupercellModel.git
cd SupercellModel
make
make vulkan
```

### Run Simulation

```bash
# Fast smoke run
./bin/tornado_sim --headless --config=configs/classic.yaml --duration=60

# Longer example run
./bin/tornado_sim --config=configs/classic.yaml \
                  --duration=3600 \
                  --write-every=30 \
                  --outdir=data/supercell_run
```

### Run Vulkan Viewer

```bash
# Vulkan bootstrap smoke check
./bin/vulkan_viewer --dry-run

# Windowed volume render
./bin/vulkan_viewer --window-test --render-backend volume --input data/supercell_run --field theta
```

### Validation Commands

```bash
make test
make test-guards
make test-backend-physics
make test-soundings
make test-radiation-regression
make test-terrain-regression
```

---

## Contributing

SupercellModel welcomes contributions from atmospheric scientists and computational researchers. The codebase is designed for extensibility and systematic comparison of parameterizations.

### Getting Started
1. **Fork** the repository on GitHub
2. **Clone** your fork: `git clone https://github.com/yourusername/SupercellModel.git`
3. **Create** a feature branch: `git checkout -b feature/your-enhancement`
4. **Build** and test: `make && ./bin/tornado_sim --config=configs/classic.yaml --duration=60`

### Development Guidelines
- **Code Style**: Modern C++17 practices, RAII, smart pointers
- **Testing**: Add validation tests for new physics parameterizations
- **Documentation**: Update relevant READMEs and `foundationalScience.md`
- **Physics**: Include scientific references for new parameterizations

### Pull Request Process
1. **Implement** your changes with comprehensive tests
2. **Validate** against existing test cases
3. **Document** scientific foundations and implementation details
4. **Submit** PR with clear description of changes and motivation

### Areas for Contribution
- **New Physics Schemes**: Additional microphysics, PBL, or turbulence options
- **Radar Operators**: Enhanced polarimetric variables or advanced sampling
- **Terrain Effects**: Terrain-science benchmarking and validation against reference cases
- **Stochastic Physics**: Ensemble workflow expansion and perturbation calibration
- **Performance**: GPU acceleration, parallelization improvements
- **Visualization**: Enhanced rendering features or new field diagnostics

---

## License & Attribution

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Scientific Foundation
SupercellModel builds upon decades of atmospheric modeling research. Key foundational references:

**Storm Dynamics & Modeling:**
- Klemp & Wilhelmson (1978) - Compressible storm dynamics
- Weisman & Klemp (1982) - Supercell simulation foundations
- Bryan et al. (2003) - Model resolution requirements

**Physics Parameterizations:**
- Kessler (1969) - Warm-rain microphysics
- Thompson et al. (2008) - Aerosol-aware microphysics
- Hong et al. (2006) - YSU boundary layer scheme

**Radar & Observations:**
- Snyder & Zhang (2003) - Radar data assimilation
- Jung et al. (2008) - Radar forward operators

For complete scientific attribution, see **[foundationalScience.md](docs/foundationalScience.md)**.

---

## Acknowledgments

SupercellModel represents a synthesis of modern atmospheric modeling techniques with an emphasis on modularity, validation, and visualization. The framework is designed to serve as both a research tool and an educational resource for understanding severe convective storms.

*Built with modern C++17 and a Vulkan-first rendering path. Inspired by the atmospheric modeling community's commitment to reproducible, well-validated research.*
