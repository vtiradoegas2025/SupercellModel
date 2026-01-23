#!/usr/bin/env python3
"""
Complete 3D supercell visualization pipeline runner.
This script demonstrates the full workflow from data conversion to video rendering.
"""
import subprocess
import sys
from pathlib import Path
import argparse
import os

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n=== {description} ===")
    print(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print("Success")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed: {e}")
        print(f"Error output: {e.stderr}")
        return False

def check_dependencies():
    """Check if required dependencies are available"""
    print("=== Checking Dependencies ===")

    required_modules = ['numpy', 'zarr', 'moderngl', 'moderngl_window', 'tqdm', 'PIL']

    missing = []

    for module in required_modules:
        try:
            if module == 'PIL':
                import PIL
            else:
                __import__(module)
            print(f"Found {module}")
        except ImportError:
            print(f"Missing {module}")
            missing.append(module)

    if missing:
        print(f"\nMissing dependencies: {', '.join(missing)}")
        print("Install with: pip install " + " ".join(missing))
        return False

    return True

def run_pipeline(args):
    """Run the complete visualization pipeline"""

    # Check dependencies first
    if not check_dependencies():
        print("Please install missing dependencies and try again.")
        return False

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data"

    success = True

    # Step 1: Convert simulation data
    if args.convert_data:
        print("\n=== Step 1: Data Conversion ===")
        if (data_dir / "exports").exists():
            cmd = f"cd {base_dir} && python visualization/convert_npy_to_zarr.py --input data/exports --output data/supercell.zarr --max_timesteps {args.max_timesteps}"
            success &= run_command(cmd, "Converting NPY data to Zarr format")
        else:
            print("No export data found. Run simulation first:")
            print("  ./bin/tornado_sim --config configs/low_res_test.yaml --output data/exports --duration 30")

    # Step 2: Run interactive viewer
    if args.run_viewer:
        if (data_dir / "supercell.zarr").exists():
            print("\n=== Step 2: Interactive 3D Viewer ===")
            print("Launching 3D supercell viewer...")
            print("Controls: SPACE (play/pause), R (reset), W (wireframe), ↑↓ (opacity), ←→ (brightness)")
            cmd = f"cd {base_dir} && python visualization/supercell_renderer.py --input data/supercell.zarr --field {args.field}"
            try:
                subprocess.run(cmd, shell=True, check=False)
            except KeyboardInterrupt:
                pass
        else:
            print("\n=== Step 2: Interactive 3D Viewer ===")
            print("No Zarr data found. Convert data first with --convert-data")

    # Step 3: Offline rendering
    if args.render_video:
        zarr_path = data_dir / "supercell.zarr"
        if not zarr_path.exists():
            print("\n=== Step 3: Offline Video Rendering ===")
            print("No Zarr data found. Convert data first with --convert-data")
            success = False
        else:
            print("\n=== Step 3: Offline Video Rendering ===")
            output_file = args.output or f"supercell_{args.field}.mp4"
            cmd = f"cd {base_dir} && python visualization/render_supercell_3d.py --input {zarr_path} --output {output_file} --fps {args.fps} --speed {args.speed}"
            if args.duration:
                cmd += f" --duration {args.duration}"

            success &= run_command(cmd, f"Rendering video to {output_file}")

    return success

def print_header(title):
    """Print a formatted header"""
    print(f"\n{'='*60}")
    print(f" {title}")
    print(f"{'='*60}")

def check_file_exists(filepath, description):
    """Check if a file exists and provide guidance"""
    if filepath.exists():
        print(f"Found {description}: {filepath}")
        return True
    else:
        print(f"Missing {description}: {filepath}")
        return False

def check_dependencies():
    """Check if required dependencies are available"""
    print("Checking Dependencies")

    required_modules = ['numpy', 'zarr', 'moderngl', 'moderngl_window', 'tqdm', 'PIL']

    missing = []

    for module in required_modules:
        try:
            if module == 'PIL':
                import PIL
            else:
                __import__(module)
            print(f"Found {module}")
        except ImportError:
            print(f"Missing {module}")
            missing.append(module)

    if missing:
        print(f"\nMissing dependencies: {', '.join(missing)}")
        print("Install with: pip install " + " ".join(missing))
        return False

    return True

def run_quick_start():
    """Run the interactive quick start guide"""
    print_header("TornadoModel 3D Visualization - Quick Start")

    print("""
Welcome to TornadoModel's 3D visualization system!

This quick start guide will help you:
1. Check if your simulation data is ready
2. Convert data to visualization format (if needed)
3. Launch an interactive 3D viewer
4. Create a sample video

Let's get started!
    """)

    # Check current directory
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    print(f"Working directory: {project_root}")

    print_header("Step 1: Check Your Simulation Data")

    # Check for simulation outputs
    exports_dir = Path("data/exports")
    if exports_dir.exists():
        npy_files = list(exports_dir.glob("**/step_*/**/*.npy"))
        if npy_files:
            print(f"Found {len(npy_files)} NPY files in {exports_dir}")
            print("   Your simulation data is ready for visualization!")

            # Show a few example files
            for i, npy_file in enumerate(npy_files[:5]):
                print(f"   Example: {npy_file.relative_to(project_root)}")
            if len(npy_files) > 5:
                print(f"   ... and {len(npy_files) - 5} more files")
        else:
            print("Found exports directory but no NPY files")
    else:
        print("No simulation data found in data/exports/")
        print("   Run a simulation first: ./bin/tornado_sim --config configs/classic.yaml --write-every 5")

    print_header("Step 2: Check Visualization Dependencies")

    if not check_dependencies():
        print("\nAfter installing dependencies, re-run this script.")
        return

    print_header("Step 3: Convert Data to Visualization Format")

    # Check if Zarr data already exists
    zarr_path = Path("data/supercell.zarr")
    if zarr_path.exists():
        print("Zarr visualization data already exists")
        print("   Skipping conversion step")
    else:
        print("Converting NPY data to Zarr format...")
        success = run_command(
            "python visualization/convert_npy_to_zarr.py --input data/exports --output data/supercell.zarr --max-timesteps 50",
            "Converting simulation data to Zarr format"
        )
        if not success:
            print("Conversion failed. Check that simulation data exists in data/exports/")
            return

    print_header("Step 4: Interactive 3D Exploration")

    print("""
Great! Your data is ready for visualization.

The interactive 3D viewer lets you explore your supercell simulation in real-time.
Controls:
  • SPACE: Play/pause animation
  • R: Reset to beginning
  • Mouse drag: Orbit camera
  • Up/Down: Adjust opacity
  • Left/Right arrows: Adjust brightness
  • 1-4: Different visualization presets
  • W: Toggle wireframe mode
    """)

    # Ask user if they want to launch viewer
    try:
        response = input("\nLaunch interactive 3D viewer now? (y/n): ").lower().strip()
        if response in ['y', 'yes']:
            print("Launching 3D viewer...")
            print("Close the window when done to continue with video rendering.")

            # Launch viewer (this will block until user closes it)
            cmd = "python visualization/supercell_renderer.py --input data/supercell.zarr --field theta"
            try:
                subprocess.run(cmd, shell=True, check=True)
                print("Viewer closed successfully")
            except KeyboardInterrupt:
                print("Viewer interrupted by user")
            except subprocess.CalledProcessError as e:
                print(f"Viewer failed: {e}")
        else:
            print("Skipping interactive viewer")
    except (EOFError, KeyboardInterrupt):
        print("\nSkipping interactive viewer")

    print_header("Step 5: Create Sample Video")

    print("""
Now let's create a short video showcasing your supercell simulation.
This will render a 15-second video at 30 FPS.
    """)

    try:
        response = input("Create sample video now? (y/n): ").lower().strip()
        if response in ['y', 'yes']:
            success = run_command(
                "python visualization/render_supercell_3d.py --input data/supercell.zarr --output sample_video.mp4 --duration 15 --fps 30",
                "Rendering sample video"
            )
            if success:
                video_path = Path("sample_video.mp4")
                if video_path.exists():
                    size_mb = video_path.stat().st_size / (1024 * 1024)
                    print(f"Created video: {video_path} ({size_mb:.1f} MB)")
                else:
                    print("Video rendering completed (file location unknown)")
            else:
                print("Video rendering failed")
        else:
            print("Skipping video creation")
    except (EOFError, KeyboardInterrupt):
        print("\nSkipping video creation")

    print_header("Quick Start Complete!")

    print("""
Congratulations! You've successfully set up TornadoModel's 3D visualization system.

What you accomplished:
• Converted simulation data to efficient Zarr format
• Explored your supercell interactively (if chosen)
• Created a sample video (if chosen)

Next Steps:
• Run more simulations with different configurations
• Experiment with different fields: --field qr, --field qv, --field qc
• Customize visualization parameters for your research needs
• Check visualization/README.md for advanced options

For scientific background, see docs/foundationalScience.md
For detailed usage, see visualization/README.md

Happy visualizing!
    """)

def main():
    parser = argparse.ArgumentParser(
        description="Complete 3D supercell visualization pipeline - bringing atmospheric simulations to life!",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
TornadoModel 3D Visualization Pipeline
======================================

This tool provides a complete workflow for visualizing supercell thunderstorm simulations:

1. Convert raw simulation data to optimized Zarr format
2. Interactive 3D exploration with real-time volume rendering
3. High-quality video generation for presentations/publications

Quick Start Examples:
  # First time setup - convert your simulation data
  python visualization/run_3d_pipeline.py --convert-data --max-timesteps 100

  # Explore interactively (requires Zarr data)
  python visualization/run_3d_pipeline.py --run-viewer --field qr

  # Create presentation video (requires Zarr data)
  python visualization/run_3d_pipeline.py --render-video --output supercell_demo.mp4 --duration 30

  # Full automated pipeline (convert → explore → render)
  python visualization/run_3d_pipeline.py --all --output my_storm.mp4

Advanced Examples:
  # Visualize different atmospheric fields
  python visualization/run_3d_pipeline.py --run-viewer --field qv      # Water vapor
  python visualization/run_3d_pipeline.py --run-viewer --field qc      # Cloud water
  python visualization/run_3d_pipeline.py --run-viewer --field theta   # Temperature

  # High-quality rendering for publications
  python visualization/run_3d_pipeline.py --render-video --output publication.mp4 --fps 60 --duration 60

Data Requirements:
- Input: Simulation NPY exports in data/exports/step_*/ directories
- Output: Zarr format for efficient 3D visualization

For detailed documentation, see:
- visualization/README.md (usage guide)
- docs/foundationalScience.md (scientific foundation)
        """
    )

    parser.add_argument('--all', action='store_true',
                       help='Run complete pipeline')
    parser.add_argument('--quick-start', action='store_true',
                       help='Interactive guided setup for new users')
    parser.add_argument('--convert-data', action='store_true',
                       help='Convert NPY exports to Zarr')
    parser.add_argument('--run-viewer', action='store_true',
                       help='Launch interactive 3D viewer')
    parser.add_argument('--render-video', action='store_true',
                       help='Render offline video')

    parser.add_argument('--max-timesteps', type=int, default=50,
                       help='Maximum timesteps to convert')
    parser.add_argument('--field', default='theta',
                       help='Primary field to visualize (theta, qr, qv, qc, qi, qs, qh, etc.)')
    parser.add_argument('--output', help='Output video file')
    parser.add_argument('--fps', type=int, default=30,
                       help='Video frame rate')
    parser.add_argument('--speed', type=float, default=1.0,
                       help='Animation speed multiplier')
    parser.add_argument('--duration', type=float,
                       help='Video duration in seconds')

    args = parser.parse_args()

    # Handle quick start mode
    if args.quick_start:
        run_quick_start()
        sys.exit(0)

    # Handle --all flag
    if args.all:
        args.convert_data = True
        args.run_viewer = True
        args.render_video = True

    # Default to viewer if no args
    if not any([args.convert_data, args.run_viewer, args.render_video]):
        args.run_viewer = True

    print("=== 3D Supercell Visualization Pipeline ===")
    print("Bringing your atmospheric simulation to life in 3D!")

    success = run_pipeline(args)

    if success:
        print("\n=== Pipeline Complete ===")
        print("3D visualization pipeline executed successfully!")
        print("Check the visualization/README.md and docs/foundationalScience.md for detailed information.")
    else:
        print("\n=== Pipeline Failed ===")
        print("Some steps failed. Check output above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()
