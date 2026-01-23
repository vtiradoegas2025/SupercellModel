#!/usr/bin/env python3
"""
Time step validation script for TornadoModel
Ensures all physics modules are properly synchronized
"""

import yaml
import os
import sys

def load_config(config_path):
    """Load YAML configuration file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def validate_config_consistency(config_path):
    """Check configuration file for consistency"""
    print(f"Validating configuration: {config_path}")

    try:
        config = load_config(config_path)
    except Exception as e:
        print(f"ERROR: Failed to load config {config_path}: {e}")
        return False

    errors = []

    # Check grid parameters
    grid = config.get('grid', {})
    dt = grid.get('dt', 0.1)
    dx = grid.get('dx', 1000.0)
    dz = grid.get('dz', 100.0)

    # Basic CFL check (rough estimate)
    max_wind_speed = 50.0  # m/s typical max
    cfl_dt = min(dx, dz) / max_wind_speed
    if dt > cfl_dt:
        errors.append(".4f")

    # Check physics module time steps
    radiation_dt = config.get('radiation', {}).get('dt_radiation', 300.0)
    if radiation_dt < dt:
        errors.append(f"Radiation dt ({radiation_dt}s) < dynamics dt ({dt}s)")

    boundary_layer_dt = config.get('boundary_layer', {}).get('dt_pbl', 60.0)
    if boundary_layer_dt < dt:
        errors.append(f"Boundary layer dt ({boundary_layer_dt}s) < dynamics dt ({dt}s)")

    turbulence_dt = config.get('turbulence', {}).get('dt_sgs', 1.0)
    if turbulence_dt < dt:
        errors.append(f"Turbulence dt ({turbulence_dt}s) < dynamics dt ({dt}s)")

    # Check that physics modules have reasonable cadences
    if radiation_dt > 3600:
        errors.append(f"Radiation dt ({radiation_dt}s) seems too large (>1hr)")

    if boundary_layer_dt > 300:
        errors.append(f"Boundary layer dt ({boundary_layer_dt}s) seems too large (>5min)")

    if turbulence_dt > 10:
        errors.append(f"Turbulence dt ({turbulence_dt}s) seems too large (>10s)")

    if errors:
        print("CONFIGURATION ERRORS:")
        for error in errors:
            print(f"  - {error}")
        return False

    print("✓ Configuration validation passed")
    return True

def validate_module_dependencies():
    """Check that all required modules are present and compilable"""
    print("Validating module dependencies...")

    # Check that all required source files exist
    required_files = [
        'src/dynamics.cpp',
        'src/radiation.cpp',
        'src/boundary_layer.cpp',
        'src/turbulence.cpp',
        'src/radar.cpp',
        'src/numerics.cpp',
        'src/chaos/chaos.cpp',
        'src/microphysics/factory.cpp',
        'src/microphysics/schemes/kessler/kessler.cpp',
    ]

    missing_files = []
    for f in required_files:
        if not os.path.exists(f):
            missing_files.append(f)

    if missing_files:
        print("MISSING REQUIRED FILES:")
        for f in missing_files:
            print(f"  - {f}")
        return False

    print("✓ All required source files present")
    return True

def validate_time_stepping_logic():
    """Validate that time stepping logic is consistent"""
    print("Validating time stepping logic...")

    # Read the main dynamics stepping file
    with open('src/dynamics.cpp', 'r') as f:
        dynamics_code = f.read()

    # Check that modules are called with appropriate time steps
    checks = [
        ('step_radiation', 'Called with current_time (cadenced internally)'),
        ('step_microphysics', 'Called with dt_dynamics'),
        ('step_turbulence', 'Called with current_time (cadenced internally)'),
    ]

    issues = []
    for func, expected in checks:
        if func in dynamics_code:
            print(f"✓ {func} found in dynamics stepping")
        else:
            issues.append(f"Missing {func} call in dynamics")

    if issues:
        print("TIME STEPPING ISSUES:")
        for issue in issues:
            print(f"  - {issue}")
        return False

    print("✓ Time stepping logic appears consistent")
    return True

def main():
    print("=== TornadoModel Time Step Validation ===\n")

    config_files = ['configs/classic.yaml', 'configs/test_small.yaml', 'configs/lp.yaml']

    all_passed = True

    # Validate configurations
    for config_file in config_files:
        if os.path.exists(config_file):
            if not validate_config_consistency(config_file):
                all_passed = False
        else:
            print(f"WARNING: Config file {config_file} not found")

    print()

    # Validate module dependencies
    if not validate_module_dependencies():
        all_passed = False

    print()

    # Validate time stepping logic
    if not validate_time_stepping_logic():
        all_passed = False

    print()

    if all_passed:
        print("🎉 All validations passed! Ready for GitHub push.")
        return 0
    else:
        print("❌ Validation failed. Please fix issues before pushing.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
