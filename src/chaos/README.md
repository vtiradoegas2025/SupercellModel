# Stochastic Perturbations (Chaos) Module

This module implements controlled stochastic variability for ensemble forecasting and uncertainty quantification in atmospheric simulations.

## Overview

Stochastic parameterizations account for:
- **Model uncertainty**: Sub-grid scale processes not explicitly resolved
- **Initial condition uncertainty**: Errors in analysis/forecasts
- **Parameter uncertainty**: Unknown physical parameters
- **Ensemble spread**: Sampling of possible weather outcomes

## Architecture

```
src/chaos/
├── chaos.cpp                # Module coordinator
├── factory.cpp/.hpp         # Chaos scheme factory
├── base/                    # Common utilities
│   ├── random_generator.cpp/.hpp    # Reproducible random numbers
│   ├── perturbation_field.cpp/.hpp  # Field perturbation utilities
│   └── correlation_filter.cpp/.hpp  # Spatial correlation
└── schemes/                 # Perturbation implementations
    ├── none/                # Deterministic (no perturbations)
    ├── initial_conditions/  # IC perturbations only
    ├── boundary_layer/      # PBL tendency perturbations
    └── full_stochastic/     # Complete SPPT implementation
```

## Stochastic Physics Fundamentals

### Stochastic Parameterized Perturbations Tendencies (SPPT)

**Concept:** Multiply parameterized tendencies by random field:

```cpp
X' = X + SPPT_factor × (X_parameterized - X_deterministic)
```

Where SPPT_factor ∈ [0,1] is a spatially/temporally correlated random field.

### Initial Condition Perturbations

**Breeding method:**
```cpp
x_{n+1} = x_deterministic + ε (x_n - x_deterministic)
```

**Ensemble transform:**
```cpp
x_ensemble = x_mean + A × x_random
```

Where A is the ensemble perturbation matrix.

### Boundary Layer Stochastic Perturbations

**Tendency perturbations:**
```cpp
∂θ/∂t|_perturbed = ∂θ/∂t|_deterministic × (1 + σ × ξ)
```

Where ξ is spatio-temporally correlated noise.

## Scheme Implementations

### None (Deterministic)

**Configuration:**
```yaml
chaos:
  scheme: "none"
```

**Features:**
- No stochastic perturbations
- Reproducible deterministic results
- Baseline for comparison studies

### Initial Conditions

**Configuration:**
```yaml
chaos:
  scheme: "initial_conditions"
  perturbation_magnitude: 0.01   # Relative perturbation amplitude
  correlation_length_km: 500.0   # Spatial correlation scale
  correlation_time_h: 24.0       # Temporal correlation
```

**Features:**
- Perturbations to initial thermodynamic fields
- Spatially correlated random fields
- Controlled ensemble spread

### Boundary Layer Tendencies

**Configuration:**
```yaml
chaos:
  scheme: "boundary_layer"
  sppt_amplitude: 0.5           # SPPT perturbation magnitude
  sppt_timescale_h: 6.0         # SPPT correlation time
  sppt_lengthscale_km: 300.0    # SPPT correlation length
  sppt_timestep_h: 1.0          # Perturbation update frequency
```

**Features:**
- Perturbations to PBL tendencies
- Time-correlated random fields
- Focus on boundary layer uncertainty

### Full Stochastic Parameterization

**Configuration:**
```yaml
chaos:
  scheme: "full_stochastic"
  sppt_amplitude: 0.3           # Overall SPPT amplitude
  sppt_timescale_h: 12.0        # Correlation time scale
  sppt_lengthscale_km: 500.0    # Correlation length scale
  stochastic_kinetics: true     # Include kinetic energy backscatter
  skebs_amplitude: 0.1          # SKEBS perturbation amplitude
```

**Features:**
- Complete SPPT implementation
- Stochastic kinetic energy backscatter (SKEBS)
- Multi-physics perturbations
- Advanced correlation modeling

## Implementation Details

### Random Number Generation

**Reproducible generation:**
```cpp
// Initialize with seed
std::mt19937 generator(seed);

// Generate Gaussian noise
std::normal_distribution<double> gaussian(0.0, 1.0);
double noise = gaussian(generator);
```

**Ensuring reproducibility:**
- Fixed seeds for ensemble members
- Deterministic sequence generation
- Platform-independent algorithms

### Spatial Correlation

**Generation of correlated fields:**
```cpp
// Fourier-based correlation
fft_noise = fft(white_noise);
correlation_filter = exp(-k² * L² / 2);
fft_correlated = fft_noise * correlation_filter;
correlated_field = ifft(fft_correlated);
```

**Correlation functions:**
- **Gaussian**: exp(-r²/(2L²))
- **Exponential**: exp(-r/L)
- **Power-law**: 1/(1 + r/L)^α

### Temporal Correlation

**AR(1) process:**
```cpp
ξ_{n+1} = α ξ_n + √(1-α²) ε_n
```

Where:
- α = exp(-Δt/τ): Autocorrelation coefficient
- τ: Correlation time scale
- ε_n: White noise

### Perturbation Application

**SPPT multiplication:**
```cpp
// Apply to physics tendencies
tendency_perturbed = tendency_deterministic * (1 + sppt_factor * noise)
```

**SKEBS addition:**
```cpp
// Add to momentum equations
du/dt += skebs_amplitude * ∇²(kinetic_energy * noise)
```

## Configuration

### Basic Ensemble Setup
```yaml
chaos:
  scheme: "initial_conditions"
  seed: 12345                 # Random seed for reproducibility
  n_ensemble_members: 20      # Number of ensemble members
```

### Advanced Stochastic Options
```yaml
chaos:
  # Spatial correlation
  correlation_model: "gaussian"  # gaussian, exponential, power_law
  horizontal_correlation_km: 300.0
  vertical_correlation_km: 2.0

  # Temporal correlation
  temporal_correlation_h: 6.0
  update_frequency_h: 1.0

  # Physical constraints
  clip_perturbations: true    # Prevent extreme values
  max_relative_change: 0.5    # Maximum 50% change
```

## Validation and Testing

### Ensemble Statistics
```bash
# Test ensemble spread
python tests/test_ensemble.py --scheme sppt --n_members 20

# Validate correlation scales
python tests/test_correlation.py --length_scale 500km --time_scale 6h
```

### Benchmark Cases
- **Perfect model**: Known truth, quantify uncertainty
- **Real cases**: Compare with operational ensembles
- **Sensitivity**: Test different perturbation amplitudes

### Expected Results

#### Ensemble Spread
- **IC perturbations**: Initial spread grows with time
- **SPPT**: Continuous spread generation and maintenance
- **SKEBS**: Enhanced spread in convective regions

#### Forecast Improvement
- **Probabilistic skill**: Better calibrated probabilities
- **Spread-skill relationship**: Spread proportional to uncertainty
- **Extreme events**: Better prediction of rare outcomes

## Performance Characteristics

### Computational Cost
- **None**: No overhead (baseline)
- **IC only**: O(N_grid) initialization
- **SPPT**: O(N_grid) per physics time step
- **Full stochastic**: O(N_grid × N_physics) per time step

### Memory Requirements
- **Random fields**: 3D arrays for spatial perturbations
- **Correlation data**: Additional arrays for filtering
- **Ensemble state**: Multiple copies for ensemble members

### Parallel Considerations
- **Ensemble members**: Independent simulations
- **Random seeds**: Unique per member
- **Load balancing**: Identical computational cost per member

## Scientific References

### Stochastic Parameterizations
- Palmer (2001) - A nonlinear dynamical perspective on model error
- Buizza et al. (1999) - Stochastic representation of model uncertainties
- Shutts (2005) - A kinetic energy backscatter algorithm

### SPPT Development
- Buizza et al. (1999) - Stochastic simulation of model uncertainties
- Palmer et al. (2009) - Stochastic parametrization and model uncertainty
- Tennant et al. (2011) - Ensemble prediction with SPPT

### Ensemble Methods
- Leith (1974) - Theoretical skill of Monte Carlo forecasts
- Epstein (1969) - Stochastic dynamic prediction
- Molteni et al. (1996) - ECMWF ensemble prediction system

This module enables probabilistic forecasting and uncertainty quantification for atmospheric simulations.
