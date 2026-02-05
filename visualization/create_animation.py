#!/usr/bin/env python3
"""
Create a simple 2D animation from NPY simulation data.
Works directly with NPY files without requiring zarr/modernGL.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import glob
import argparse
from tqdm import tqdm
try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False
    print("Warning: PIL not available")
try:
    import imageio
    IMAGEIO_AVAILABLE = True
except ImportError:
    IMAGEIO_AVAILABLE = False

def load_timestep_data(step_dir, field_name, NR=256, NTH=128, NZ=128):
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
                return None
        else:
            return None
    
    if len(data_slices) != NTH:
        return None
    
    # Stack into 3D array (NZ, NTH, NR)
    # Each slice is (NZ, NR), stacking along axis=1 gives (NZ, NTH, NR)
    stacked = np.stack(data_slices, axis=1)
    # Transpose to (NR, NTH, NZ) for easier indexing
    return np.transpose(stacked, (2, 1, 0))

def create_horizontal_slice(data_3d, z_level=None):
    """
    Create a horizontal slice from 3D cylindrical data.
    If z_level is None, uses middle level.
    
    Args:
        data_3d: 3D array (NR, NTH, NZ)
        z_level: Vertical level index (None = middle)
    
    Returns:
        2D array (NR, NTH) - horizontal slice
    """
    if z_level is None:
        z_level = data_3d.shape[2] // 2
    
    return data_3d[:, :, z_level]

def cylindrical_to_cartesian_2d(data_cylindrical, NR, NTH):
    """
    Convert 2D cylindrical slice to Cartesian coordinates.
    
    Args:
        data_cylindrical: 2D array (NR, NTH)
        NR, NTH: Grid dimensions
    
    Returns:
        data_cartesian: 2D array (NX, NY)
        x_coords, y_coords: coordinate arrays
    """
    # Create radial and azimuthal coordinates
    dr = 1000.0  # meters
    r_coords = np.arange(NR) * dr
    theta_coords = np.linspace(0, 2*np.pi, NTH, endpoint=False)
    
    # Target Cartesian grid
    NX, NY = 128, 128
    x_max = r_coords[-1]
    y_max = r_coords[-1]
    
    x_coords = np.linspace(-x_max, x_max, NX)
    y_coords = np.linspace(-y_max, y_max, NY)
    
    # Create meshgrid
    X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')
    
    # Convert to cylindrical
    R = np.sqrt(X**2 + Y**2)
    Theta = np.arctan2(Y, X)
    Theta = np.where(Theta < 0, Theta + 2*np.pi, Theta)
    
    # Interpolate
    data_cartesian = np.zeros((NX, NY))
    
    r_indices = np.searchsorted(r_coords, R.ravel()) - 1
    theta_indices = np.searchsorted(theta_coords, Theta.ravel()) - 1
    
    r_indices = np.clip(r_indices, 0, NR-1)
    theta_indices = np.clip(theta_indices, 0, NTH-1)
    
    data_cartesian.ravel()[:] = data_cylindrical[r_indices, theta_indices]
    
    return data_cartesian, x_coords, y_coords

def main():
    parser = argparse.ArgumentParser(description="Create animation from NPY simulation data")
    parser.add_argument("--input", required=True, help="Input directory with step_*/ subdirs")
    parser.add_argument("--output", default="animation.mp4", help="Output video file")
    parser.add_argument("--field", default="theta", help="Field to visualize (theta, qr, qv, etc.)")
    parser.add_argument("--z-level", type=int, help="Vertical level (default: middle)")
    parser.add_argument("--fps", type=int, default=10, help="Frames per second")
    parser.add_argument("--max-timesteps", type=int, help="Maximum timesteps to process")
    parser.add_argument("--dpi", type=int, default=100, help="DPI for output")
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    if not input_dir.exists():
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    # Find all timestep directories
    step_dirs = sorted(glob.glob(str(input_dir / "step_*")))
    if not step_dirs:
        raise ValueError(f"No step directories found in {input_dir}")
    
    max_timesteps = args.max_timesteps or len(step_dirs)
    step_dirs = step_dirs[:max_timesteps]
    
    print(f"Found {len(step_dirs)} timesteps")
    print(f"Field: {args.field}, Z-level: {args.z_level or 'middle'}")
    
    # Grid parameters (matching config)
    NR, NTH, NZ = 256, 128, 128
    
    # Load all timesteps
    print("Loading data...")
    timestep_data = []
    for step_dir in tqdm(step_dirs, desc="Loading"):
        data_3d = load_timestep_data(step_dir, args.field, NR, NTH, NZ)
        if data_3d is None:
            print(f"Warning: Failed to load {args.field} from {step_dir}")
            continue
        
        # Extract horizontal slice
        slice_2d = create_horizontal_slice(data_3d, args.z_level)
        
        # Convert to Cartesian
        data_cart, x_coords, y_coords = cylindrical_to_cartesian_2d(slice_2d, NR, NTH)
        timestep_data.append(data_cart)
    
    if not timestep_data:
        raise ValueError("No data loaded!")
    
    print(f"Loaded {len(timestep_data)} timesteps")
    
    # Create animation by saving frames and combining
    print("Creating animation frames...")
    
    # Determine color scale
    all_data = np.concatenate([d.ravel() for d in timestep_data])
    vmin = np.nanpercentile(all_data, 1)
    vmax = np.nanpercentile(all_data, 99)
    
    # Create temporary directory for frames
    temp_dir = Path("/tmp") / f"animation_frames_{Path(args.output).stem}"
    temp_dir.mkdir(exist_ok=True)
    
    # Save individual frames
    frame_files = []
    plt.ioff()
    for i, data in enumerate(tqdm(timestep_data, desc="Creating frames")):
        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(data, 
                       extent=[x_coords[0]/1000, x_coords[-1]/1000, 
                               y_coords[0]/1000, y_coords[-1]/1000],
                       origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_xlabel('X (km)', fontsize=12)
        ax.set_ylabel('Y (km)', fontsize=12)
        ax.set_title(f'{args.field} at z={args.z_level or NZ//2} (Timestep {i})', fontsize=14)
        plt.colorbar(im, ax=ax, label=args.field)
        
        frame_file = temp_dir / f"frame_{i:04d}.png"
        plt.savefig(frame_file, dpi=args.dpi, bbox_inches='tight')
        plt.close(fig)
        frame_files.append(str(frame_file))
    
    # Combine frames into animation
    if args.output.endswith('.gif') and PIL_AVAILABLE:
        print(f"Combining {len(frame_files)} frames into GIF...")
        try:
            images = [Image.open(f) for f in frame_files]
            images[0].save(args.output, save_all=True, append_images=images[1:], 
                          duration=1000/args.fps, loop=0)
            print(f"Animation saved to {args.output}")
            # Cleanup
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
        except Exception as e:
            print(f"Error creating GIF: {e}")
            print(f"Frames saved in {temp_dir}")
    elif IMAGEIO_AVAILABLE:
        print(f"Combining {len(frame_files)} frames into {args.output}...")
        try:
            if args.output.endswith('.gif'):
                imageio.mimsave(args.output, [imageio.imread(f) for f in frame_files], 
                               fps=args.fps, loop=0)
            else:
                imageio.mimsave(args.output, [imageio.imread(f) for f in frame_files], 
                               fps=args.fps, codec='libx264', quality=8)
            print(f"Animation saved to {args.output}")
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
        except Exception as e:
            print(f"Error combining frames: {e}")
            print(f"Frames saved in {temp_dir}")
    else:
        print(f"Animation libraries not available. Frames saved in {temp_dir}")
        print(f"To create animation:")
        print(f"  For GIF: Use PIL or imageio")
        print(f"  For MP4: ffmpeg -r {args.fps} -i {temp_dir}/frame_%04d.png -vcodec libx264 {args.output}")

if __name__ == "__main__":
    main()
