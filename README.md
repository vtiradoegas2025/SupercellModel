# SupercellModel

[![C++](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](#license)

A high-performance atmospheric simulation framework for **supercell and tornadic storm research**, featuring a modular physics design (**C++17**) and a companion **3D visualization pipeline** (**Python + OpenGL**).

> **Status:** Work in progress / research prototype. APIs, file formats, and results may change.

---

## Quick Start

```bash
# Clone and build
git clone https://github.com/vtiradoegas2025/SupercellModel.git
cd SupercellModel && make

# Run a test simulation
./bin/tornado_sim --headless --config=configs/classic.yaml --duration=60

# Visualize results
pip install -r requirements.txt
python visualization/run_3d_pipeline.py --quick-start
```

---

## What is SupercellModel?

SupercellModel is a research-grade atmospheric simulation framework designed for studying severe convective storms. It implements the compressible, non-hydrostatic equations in cylindrical coordinates with a modular physics architecture that enables systematic comparison of parameterizations.

### Key Features

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
- **Radiation**: 2 schemes (Simple grey, RRTMG)
- **Radar Operators**: Forward simulators for Z/Ze, V_r, Z_DR, K_DP, ρ_HV

#### 3D Visualization Pipeline (Python + OpenGL)
- **Volume Rendering**: GPU-accelerated ray marching through atmospheric volumes
- **Coordinate Transformation**: Cylindrical simulation data → Cartesian visualization
- **Interactive Viewer**: Real-time exploration with orbital camera controls
- **Offline Rendering**: High-quality video production with FFmpeg

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
│   ├── dynamics/          # Compressible Euler solver
│   ├── microphysics/      # Cloud physics parameterizations
│   ├── boundary_layer/    # Planetary boundary layer schemes
│   ├── turbulence/        # Sub-grid turbulence closures
│   ├── radiation/         # Radiative transfer schemes
│   ├── radar/            # Forward radar operators
│   └── numerics/         # Advection/diffusion/time-stepping
├── include/               # Public headers and API
├── configs/               # Pre-configured simulation setups
├── visualization/         # 3D rendering pipeline (Python)
├── docs/                  # Technical documentation
├── tests/                 # Validation and integration tests
├── bin/                   # Build outputs
└── data/                  # Sample data and exports
```

---

## Prerequisites & Installation

### System Requirements
- **C++ Compiler**: C++17 compliant (clang++ ≥ 9.0, g++ ≥ 7.0)
- **OpenMP**: Optional but recommended for multi-core performance (libomp on macOS via Homebrew)
- **Python**: 3.10+ with scientific stack
- **Graphics**: OpenGL 3.3+ compatible GPU (for visualization)
- **Memory**: 8GB+ RAM recommended for production simulations
- **CPU**: Multi-core processor recommended (4+ cores for optimal performance)

### Dependencies
- **C++**: SFML 3.x (optional GUI), Standard Template Library, OpenMP (optional)
- **Python**: numpy, zarr, moderngl, pillow, tqdm, ffmpeg

### Build & Install

```bash
# Clone repository
git clone https://github.com/vtiradoegas2025/SupercellModel.git
cd SupercellModel

# Build simulation engine (with optimizations enabled)
make

# Install Python dependencies
pip install -r requirements.txt

# Optional: Build with GUI support
make GUI=1

# Note: OpenMP is automatically detected and enabled if available
# On macOS, install libomp via: brew install libomp
```

---

## Usage Examples

### Basic Simulation Run

```bash
# Run classic supercell case
./bin/tornado_sim --config=configs/classic.yaml \
                  --duration=3600 \
                  --write-every=30 \
                  --outdir=data/supercell_run

# Config file specifies:
# - CAPE = 2500 J/kg, shear = 40 m/s
# - Thompson microphysics, YSU PBL, Smagorinsky turbulence
# - Domain: 256×128×128 grid (Δr = 500m, Δz = 100m)
```

### Interactive Visualization

```bash
# Launch 3D viewer
python visualization/supercell_renderer.py \
    --input data/supercell_run \
    --field theta

# Controls:
# SPACE: play/pause, R: reset, W: wireframe toggle
# Mouse: orbit camera, Up/Down: opacity, Left/Right: brightness
```

### Automated Pipeline

```bash
# Complete simulation + visualization workflow
python visualization/run_3d_pipeline.py --all \
    --config configs/classic.yaml \
    --output supercell_showcase.mp4
```

### Configuration Examples

**Classic Supercell** (Weisman-Klemp):
```yaml
microphysics: thompson
boundary_layer: ysu
turbulence: smagorinsky
cape_jkg: 2500
shear_ms: 40
domain_km: 128
```

**High-Resolution Tornado Case**:
```yaml
grid:
  nr: 512
  nth: 256
  nz: 256
microphysics: milbrandt
radar: full_suite
```

---

## Physics Modules

### Microphysics Parameterizations
| Scheme | Species | Features | Use Case |
|--------|---------|----------|----------|
| **Kessler** | q_v, q_c, q_r | Warm-rain, autoconversion | Idealized studies |
| **Thompson** | q_v, q_c, q_r, q_i, q_s, q_g | Aerosol-aware, ice processes | Realistic storms |
| **Lin** | q_v, q_c, q_r, q_i, q_s | Fixed ice properties | Mixed-phase research |
| **Milbrandt-Yau** | q_v + 6 hydrometeors | Double-moment | Size distribution studies |

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
- Terrain module currently excluded from build
- Chaos/stochastic perturbations are placeholders
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
- **[Validation Framework](tests/README.md)** - Physics and numerical checks
- **[Visualization Guide](visualization/README.md)** - 3D rendering pipeline

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
- **Terrain Effects**: Complete terrain integration and testing
- **Stochastic Physics**: Full SPPT and ensemble capability
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

*Built with modern C++17 and scientific Python. Inspired by the atmospheric modeling community's commitment to reproducible, well-validated research.*
