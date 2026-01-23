#!/usr/bin/env python3
"""
Quick Start Guide for TornadoModel 3D Visualization

This script provides a simple, guided introduction to visualizing your supercell simulations.
Run this script to get started with 3D visualization of your atmospheric data.

Usage:
    python visualization/quick_start.py
"""
import os
import sys
from pathlib import Path
import subprocess

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

def run_command(cmd, description, cwd=None):
    """Run a command and return success status"""
    print(f"\nRunning: {description}")
    print(f"   Command: {cmd}")

    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, check=True)
        print("   Success!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"   Failed: {e}")
        if e.stderr:
            print(f"   Error: {e.stderr}")
        return False

def main():
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

    # Check Python dependencies
    required_modules = ['numpy', 'zarr', 'moderngl', 'moderngl_window', 'tqdm', 'PIL']
    missing_modules = []

    for module in required_modules:
        try:
            if module == 'PIL':
                import PIL
            else:
                __import__(module)
            print(f"Found {module}")
        except ImportError:
            print(f"Missing {module}")
            missing_modules.append(module)

    if missing_modules:
        print(f"\nMissing dependencies: {', '.join(missing_modules)}")
        print("Install with: pip install " + " ".join(missing_modules))
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
  • ↑↓: Adjust opacity
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

if __name__ == "__main__":
    main()
