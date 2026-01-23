# Atmospheric Soundings Module

This module provides support for loading and interpolating atmospheric sounding data from external sources, enabling the use of real atmospheric profiles for model initialization instead of procedural profiles.

## Architecture

The soundings module follows the same factory pattern as other SupercellModel components (microphysics, terrain, radiation, etc.):

```
src/soundings/
├── base/                          # Base classes and utilities
│   ├── soundings_base.hpp/cpp     # Data structures and base class
├── schemes/                       # Sounding implementations
│   └── sharpy/                    # SHARPY sounding reader
│       ├── sharpy_sounding.hpp/cpp
├── factory.hpp/cpp                # Factory pattern implementation
├── soundings.cpp                  # Main API
├── integration_example.hpp        # Integration guide
└── README.md                      # This file
```

## Features

- **Multiple data formats**: Support for SHARPY HDF5 and NetCDF files
- **Quality control**: Automatic validation and filtering of sounding data
- **Interpolation**: Linear interpolation to model grid heights
- **Fallback support**: Graceful degradation to procedural profiles
- **Extensible design**: Easy to add new sounding formats

## Usage

### Basic Usage

```cpp
#include "soundings.hpp"

// Configure sounding system
SoundingConfig config;
config.scheme_id = "sharpy";
config.file_path = "path/to/sounding.h5";
config.use_fallback_profiles = true;

// Initialize
initialize_soundings(config);

// Load sounding data
SoundingData data = load_sounding_data();

// Interpolate to model grid
std::vector<double> model_heights = {0, 100, 200, ..., 15000};
SoundingData interpolated = interpolate_sounding_to_grid(data, model_heights);
```

### Configuration

Add to YAML config files:

```yaml
environment:
  sounding:
    scheme_id: "sharpy"           # "sharpy", "none", or empty
    file_path: "data/sounding.h5"  # Path to sounding file
    use_fallback_profiles: true    # Use procedural if sounding fails
    interpolation_method: 0        # 0=linear, 1=spline, 2=log-linear
    extrapolate_below_ground: false
    extrapolate_above_top: false
```

## SHARPY Integration

SHARPY is a Python package for analyzing atmospheric sounding data. This module can read SHARPY-exported files containing:

- Height/pressure profiles
- Temperature and dewpoint
- Wind speed and direction
- Derived quantities (potential temperature, mixing ratio, etc.)

### Expected SHARPY HDF5 Structure

```
/profiles/
  height_m                 # Height above ground (m)
  pressure_hpa            # Pressure (hPa)
  temperature_k           # Temperature (K)
  dewpoint_k              # Dewpoint temperature (K)
  wind_speed_ms           # Wind speed (m/s)
  wind_direction_deg      # Wind direction (degrees)
/metadata/
  station_id              # Station identifier
  timestamp_utc           # Observation time
  latitude_deg            # Station latitude
  longitude_deg           # Station longitude
  elevation_m             # Station elevation
```

## Implementation Status

### Current Implementation
- ✅ Data structures and base classes
- ✅ Factory pattern and API
- ✅ Quality control and validation
- ✅ Linear interpolation
- ✅ Sample sounding generator for testing
- ✅ Integration examples and documentation

### Placeholder (Needs Library Dependencies)
- 🔄 HDF5 file reading (requires libhdf5)
- 🔄 NetCDF file reading (requires libnetcdf)
- 🔄 Spline interpolation (requires math library)

## Dependencies

### Required for Full Functionality
- **HDF5**: For reading SHARPY HDF5 files
  - Ubuntu: `sudo apt-get install libhdf5-dev`
  - macOS: `brew install hdf5`
- **NetCDF**: For reading NetCDF files
  - Ubuntu: `sudo apt-get install libnetcdf-dev`
  - macOS: `brew install netcdf`

### Build Integration

Add to Makefile:
```makefile
SRC += src/soundings/soundings.cpp \
       src/soundings/factory.cpp \
       src/soundings/base/soundings_base.cpp \
       src/soundings/schemes/sharpy/sharpy_sounding.cpp

INCLUDE += -Iinclude

# Add when libraries are available:
# LIBS += -lhdf5 -lnetcdf
```

## Testing

The module includes a sample sounding generator for testing without external files. To test with real SHARPY data:

1. Install SHARPY: `pip install sharpy`
2. Create a sounding file in Python:
   ```python
   import sharpy
   import numpy as np

   # Load or create sounding data
   sounding = sharpy.Sounding(...)  # Your sounding data

   # Export to HDF5
   sounding.to_hdf('sounding.h5')
   ```
3. Use in SupercellModel configuration

## Future Enhancements

- Additional interpolation methods (spline, log-linear)
- Support for multiple sounding sources
- Sounding data assimilation
- Real-time sounding updates
- Additional quality control checks
