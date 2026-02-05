#!/usr/bin/env python3
"""
Convert NPY export sequence to Zarr format with cylindrical to Cartesian transformation.
"""
import numpy as np
import glob
import argparse
from pathlib import Path
import json
from tqdm import tqdm

# Import zarr only when needed to avoid import errors during help
try:
    import zarr
    ZARR_AVAILABLE = True
except ImportError:
    ZARR_AVAILABLE = False

def load_grid_parameters():
    """Load grid parameters from the simulation (these should match equations.cpp)"""
    # These values should match your C++ simulation
    # Updated to match classic.yaml config: nx=256, ny=128, nz=128
    NR = 256   # radial points
    NTH = 128  # azimuthal points
    NZ = 128   # vertical points
    dr = 1000.0  # radial resolution (m) - matches config grid.dx
    dz = 100.0  # vertical resolution (m) - matches config grid.dz

    # Create coordinate arrays (cylindrical)
    r_coords = np.arange(NR) * dr
    theta_coords = np.linspace(0, 2*np.pi, NTH, endpoint=False)
    z_coords = np.arange(NZ) * dz

    return NR, NTH, NZ, r_coords, theta_coords, z_coords

def cylindrical_to_cartesian_simple(r_coords, theta_coords, z_coords, data_cylindrical):
    """
    Simple cylindrical to Cartesian conversion using nearest neighbor interpolation.

    Args:
        r_coords: 1D array of radial coordinates
        theta_coords: 1D array of azimuthal angles
        z_coords: 1D array of vertical coordinates
        data_cylindrical: 3D array (NR, NTH, NZ)

    Returns:
        data_cartesian: 3D array (NX, NY, NZ)
        x_coords, y_coords, z_coords: coordinate arrays
    """
    NR, NTH, NZ = data_cylindrical.shape

    # Target Cartesian grid
    NX, NY = 128, 128  # Start with reasonable resolution
    x_max = r_coords[-1]
    y_max = r_coords[-1]
    z_max = z_coords[-1]

    x_coords = np.linspace(-x_max, x_max, NX)
    y_coords = np.linspace(-y_max, y_max, NY)
    z_coords_cart = np.linspace(0, z_max, NZ)

    # Create Cartesian grid
    X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords_cart, indexing='ij')

    # Convert to cylindrical coordinates
    R = np.sqrt(X**2 + Y**2)
    Theta = np.arctan2(Y, X)

    # Handle negative angles
    Theta = np.where(Theta < 0, Theta + 2*np.pi, Theta)

    data_cartesian = np.zeros((NX, NY, NZ))

    # Simple nearest neighbor interpolation
    r_indices = np.searchsorted(r_coords, R.ravel()) - 1
    theta_indices = np.searchsorted(theta_coords, Theta.ravel()) - 1
    z_indices = np.searchsorted(z_coords, Z.ravel()) - 1

    # Clamp indices
    r_indices = np.clip(r_indices, 0, NR-1)
    theta_indices = np.clip(theta_indices, 0, NTH-1)
    z_indices = np.clip(z_indices, 0, NZ-1)

    # Interpolate
    data_cartesian.ravel()[:] = data_cylindrical[r_indices, theta_indices, z_indices]

    return data_cartesian, x_coords, y_coords, z_coords_cart

def load_theta_slice_data(step_dir, field_name, NR, NTH, NZ):
    """
    Load all theta slices for a given field and timestep.

    Args:
        step_dir: Directory containing NPY files for this timestep
        field_name: Name of the field (e.g., 'theta', 'qr', 'u')
        NR, NTH, NZ: Grid dimensions

    Returns:
        3D array (NR, NTH, NZ) or None if data not found
    """
    data_slices = []

    for theta_idx in range(NTH):
        npy_file = Path(step_dir) / f"th{theta_idx}_{field_name}.npy"
        if npy_file.exists():
            try:
                slice_data = np.load(npy_file)
                data_slices.append(slice_data)
            except Exception as e:
                print(f"Warning: Failed to load {npy_file}: {e}")
                continue
        else:
            print(f"Warning: Missing file {npy_file}")
            return None

    if len(data_slices) != NTH:
        print(f"Warning: Expected {NTH} theta slices, got {len(data_slices)}")
        return None

    # Stack into 3D array (NR, NTH, NZ)
    return np.stack(data_slices, axis=1)

def auto_detect_fields(input_dir, NR, NTH, NZ):
    """
    Auto-detect available fields by scanning the first timestep directory.
    
    Args:
        input_dir: Path to directory containing step_*/ subdirs
        NR, NTH, NZ: Grid dimensions
    
    Returns:
        List of detected field names
    """
    step_dirs = sorted(glob.glob(str(Path(input_dir) / "step_*")))
    if not step_dirs:
        return []
    
    # Scan first timestep directory
    first_step_dir = Path(step_dirs[0])
    detected_fields = set()
    
    # Look for files matching pattern th{theta_idx}_{field_name}.npy
    for npy_file in first_step_dir.glob("th*_*.npy"):
        # Extract field name from filename like "th0_theta.npy" -> "theta"
        parts = npy_file.stem.split("_", 1)  # Split on first underscore
        if len(parts) == 2:
            field_name = parts[1]  # Everything after "th0_"
            detected_fields.add(field_name)
    
    return sorted(list(detected_fields))

def create_zarr_metadata(store, grid_params, case_name="converted"):
    """Create Zarr metadata as per supercell spec"""
    NR, NTH, NZ, r_coords, theta_coords, z_coords = grid_params

    # Root attributes
    store.attrs['version'] = '1.0'
    store.attrs['dx'] = 100.0  # Should match your simulation
    store.attrs['dy'] = 100.0
    store.attrs['dz'] = 100.0
    store.attrs['dt'] = 0.1
    store.attrs['time_units'] = 'seconds'
    store.attrs['origin'] = [0.0, 0.0, 0.0]
    store.attrs['case_name'] = case_name

    # Coordinate arrays
    x_max = r_coords[-1]
    y_max = r_coords[-1]
    z_max = z_coords[-1]

    NX, NY = 128, 128  # Target Cartesian resolution
    x_coords = np.linspace(-x_max, x_max, NX)
    y_coords = np.linspace(-y_max, y_max, NY)
    z_coords_cart = np.linspace(0, z_max, NZ)

    store['x'] = x_coords
    store['y'] = y_coords
    store['z'] = z_coords_cart

def diagnose_zarr_structure(zarr_path):
    """Diagnose the structure of the zarr file"""
    if not ZARR_AVAILABLE:
        print("Error: zarr package not installed. Run: pip install zarr")
        return False

    zarr_file = Path(zarr_path)
    if not zarr_file.exists():
        print(f"Error: {zarr_file} not found")
        return False

    try:
        store = zarr.open(str(zarr_file), mode='r')
        print(f"=== Zarr Store Structure: {zarr_file} ===")
        print(f"Root keys: {list(store.keys())}")

        for key in store.keys():
            item = store[key]
            print(f"\n{key}:")
            print(f"  Type: {type(item)}")

            if hasattr(item, 'shape'):
                print(f"  Shape: {item.shape}")
                print(f"  Dtype: {item.dtype}")
                print(f"  Size: {item.size} elements")
            elif hasattr(item, 'keys'):
                print(f"  Keys: {list(item.keys())}")
                # Check first few sub-items
                sub_keys = list(item.keys())[:3]
                for sub_key in sub_keys:
                    sub_item = item[sub_key]
                    if hasattr(sub_item, 'shape'):
                        print(f"    {key}/{sub_key}: shape {sub_item.shape}")
                    else:
                        print(f"    {key}/{sub_key}: {type(sub_item)}")

        return True

    except Exception as e:
        print(f"Error reading zarr file: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Convert NPY exports to Zarr format")
    parser.add_argument("--input", help="Input directory with step_*/ subdirs")
    parser.add_argument("--output", help="Output Zarr file path")
    parser.add_argument("--fields", nargs="+", default=["theta", "qv", "qc", "qr", "qh", "qg", "u", "v", "w", "rho", "p", "radar", "tracer"],
                       help="Fields to convert (default includes all exported fields)")
    parser.add_argument("--auto-detect-fields", action="store_true",
                       help="Auto-detect available fields from first timestep directory")
    parser.add_argument("--validate-data", action="store_true",
                       help="Validate data quality (check for NaN, report statistics)")
    parser.add_argument("--max_timesteps", type=int, help="Maximum timesteps to convert")
    parser.add_argument("--diagnose", help="Diagnose existing Zarr file instead of converting")
    args = parser.parse_args()

    # Check if zarr is available
    if not ZARR_AVAILABLE:
        print("Error: zarr package not installed. Install with: pip install zarr")
        return

    # Handle diagnose mode
    if args.diagnose:
        if not args.output:
            print("Error: --output required for --diagnose")
            return
        diagnose_zarr_structure(args.output)
        return

    input_dir = Path(args.input)
    output_path = Path(args.output)

    if not input_dir.exists():
        raise ValueError(f"Input directory {input_dir} does not exist")

    # Get grid parameters
    NR, NTH, NZ, r_coords, theta_coords, z_coords = load_grid_parameters()

    # Find all timestep directories
    step_dirs = sorted(glob.glob(str(input_dir / "step_*")))
    if not step_dirs:
        raise ValueError(f"No step directories found in {input_dir}")

    print(f"Found {len(step_dirs)} timesteps")

    # Auto-detect fields if requested
    fields_to_convert = args.fields
    if args.auto_detect_fields:
        detected_fields = auto_detect_fields(input_dir, NR, NTH, NZ)
        if detected_fields:
            fields_to_convert = detected_fields
            print(f"Auto-detected {len(detected_fields)} fields: {', '.join(detected_fields)}")
        else:
            print("Warning: No fields detected, using default field list")
            print(f"Using default fields: {', '.join(fields_to_convert)}")

    # Create Zarr store
    store = zarr.open(str(output_path), mode='w')

    # Create metadata
    grid_params = (NR, NTH, NZ, r_coords, theta_coords, z_coords)
    create_zarr_metadata(store, grid_params, output_path.stem)

    # Convert timesteps
    max_timesteps = args.max_timesteps or len(step_dirs)
    actual_timesteps = min(max_timesteps, len(step_dirs))

    # Collect all data first
    field_data = {field: [] for field in fields_to_convert}

    for i in tqdm(range(actual_timesteps), desc="Loading data"):
        step_path = Path(step_dirs[i])
        step_dir = str(step_path)

        for field_name in fields_to_convert:
            # Load cylindrical data
            data_cylindrical = load_theta_slice_data(step_dir, field_name, NR, NTH, NZ)
            if data_cylindrical is None:
                print(f"Warning: Missing data for {field_name} at timestep {i}")
                continue

            # Convert to Cartesian
            data_cartesian, x_coords, y_coords, z_coords_cart = cylindrical_to_cartesian_simple(
                r_coords, theta_coords, z_coords, data_cylindrical
            )

            # Optional data validation
            if args.validate_data and i == 0:  # Only validate first timestep for performance
                nan_count = np.isnan(data_cartesian).sum()
                nan_fraction = nan_count / data_cartesian.size
                data_min = np.nanmin(data_cartesian)
                data_max = np.nanmax(data_cartesian)
                data_mean = np.nanmean(data_cartesian)
                print(f"  {field_name}: min={data_min:.4e}, max={data_max:.4e}, mean={data_mean:.4e}, "
                      f"NaN={nan_count} ({100*nan_fraction:.2f}%)")

            field_data[field_name].append(data_cartesian)

    # Store as single arrays
    for field_name, data_list in field_data.items():
        if data_list:
            # Stack all timesteps into a single 4D array (time, z, y, x)
            field_array = np.stack(data_list, axis=0)

            # Map field names to visualization-friendly names
            # Note: The simulation exports "radar" but visualization expects "reflectivity_dbz"
            if field_name == "radar":
                store_name = "reflectivity_dbz"
            else:
                store_name = field_name

            store[store_name] = field_array
            print(f"Stored {field_name} as {store_name}: shape {field_array.shape}")

    print(f"Conversion complete. Output: {output_path}")

if __name__ == "__main__":
    main()
