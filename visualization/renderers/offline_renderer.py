"""
Offline renderer for creating video animations from simulation data.
"""

import argparse
import subprocess
import shutil
import tempfile
from pathlib import Path
from typing import Optional
import numpy as np
from PIL import Image
from tqdm import tqdm

from .base_renderer import BaseRenderer
from ..core.transfer_function import ColorTransferFunction, GrayscaleTransferFunction
from ..core.camera import OrbitCamera, ScientificCamera
from ..passes.volume_pass import VolumeRenderPass


class OfflineRenderer(BaseRenderer):
    """
    Offline renderer for creating video animations.
    
    Renders frames to disk and encodes them into video files using FFmpeg.
    """
    
    def __init__(self, data_path: str,
                 primary_field: str = "theta",
                 shader: str = "volume_color",
                 width: int = 1920,
                 height: int = 1080):
        """
        Initialize offline renderer.
        
        Args:
            data_path: Path to simulation data
            primary_field: Primary field to visualize
            shader: Shader name (volume_color, volume_grayscale, etc.)
            width: Video width
            height: Video height
        """
        super().__init__(data_path, primary_field, width=width, height=height)
        
        # Setup volume render pass
        self.volume_pass = VolumeRenderPass(
            self.engine.shader_manager,
            self.transfer_function,
            self.volume_texture,
            vertex_shader="volume",
            fragment_shader=shader
        )
        self.engine.add_render_pass(self.volume_pass)
        
        # Determine number of timesteps
        if self.data_manager.format == 'zarr':
            # Get from data
            sample_field = self.data_manager.get_all_timesteps(self.primary_field)
            self.num_timesteps = sample_field.shape[0] if sample_field is not None else 1
        else:
            # Count step directories
            import glob
            step_dirs = glob.glob(str(self.data_path / "step_*"))
            self.num_timesteps = len(step_dirs)
    
    def render_animation(self, output_path: str,
                        fps: int = 30,
                        duration: Optional[float] = None,
                        speed: float = 1.0,
                        camera_rotation: bool = True):
        """
        Render animation to video file.
        
        Args:
            output_path: Output video file path
            fps: Frame rate
            duration: Duration in seconds (if None, render all timesteps)
            speed: Animation speed multiplier
            camera_rotation: If True, rotate camera around scene
        """
        if duration is None:
            frame_count = self.num_timesteps
        else:
            frame_count = int(duration * fps)
        
        # Create temporary directory for frames
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            print(f"Rendering {frame_count} frames...")
            
            # Render each frame
            for frame_idx in tqdm(range(frame_count), desc="Rendering frames", unit="frame"):
                timestep = int(frame_idx * speed) % self.num_timesteps
                
                # Update camera angle if rotating
                if camera_rotation and isinstance(self.camera, OrbitCamera):
                    angle = frame_idx * 0.02  # Slow rotation
                    self.camera.set_angle(angle)
                
                # Load timestep data
                self.load_timestep(timestep)
                
                # Render frame
                pixels = self.engine.render_frame(timestep=timestep)
                
                # Save as PNG
                img = Image.fromarray(pixels, 'RGBA')
                frame_path = temp_path / f"frame_{frame_idx:06d}.png"
                img.save(frame_path)
            
            print("Encoding video...")
            
            # Check if FFmpeg is available
            ffmpeg_path = shutil.which("ffmpeg")
            if ffmpeg_path is None:
                raise RuntimeError(
                    "FFmpeg is not installed or not in PATH. "
                    "Please install FFmpeg to encode videos.\n"
                    "Installation: https://ffmpeg.org/download.html"
                )
            
            # Encode to video
            cmd = [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(temp_path / "frame_%06d.png"),
                "-c:v", "libx264",
                "-preset", "slow",
                "-crf", "22",
                "-pix_fmt", "yuv420p",
                "-vf", f"scale={self.width}:{self.height}:force_original_aspect_ratio=decrease,pad={self.width}:{self.height}:(ow-iw)/2:(oh-ih)/2",
                str(output_path)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"FFmpeg error: {result.stderr}")
                raise RuntimeError("Video encoding failed")
        
        print(f"Video saved to: {output_path}")


def main():
    """Command-line interface for offline rendering"""
    parser = argparse.ArgumentParser(description="Render 3D supercell animations")
    parser.add_argument("--input", required=True, help="Zarr dataset path or NPY directory")
    parser.add_argument("--output", required=True, help="Output video file")
    parser.add_argument("--field", default="theta", help="Primary field to visualize")
    parser.add_argument("--shader", default="volume_color", 
                       choices=["volume_color", "volume_grayscale"],
                       help="Shader to use")
    parser.add_argument("--fps", type=int, default=30, help="Frame rate")
    parser.add_argument("--duration", type=float, help="Duration in seconds")
    parser.add_argument("--speed", type=float, default=1.0, help="Animation speed multiplier")
    parser.add_argument("--width", type=int, default=1920, help="Video width")
    parser.add_argument("--height", type=int, default=1080, help="Video height")
    parser.add_argument("--camera", default="orbit", choices=["orbit", "scientific"],
                       help="Camera type")
    
    args = parser.parse_args()
    
    # Create renderer
    renderer = OfflineRenderer(
        args.input,
        primary_field=args.field,
        shader=args.shader,
        width=args.width,
        height=args.height
    )
    
    # Set camera
    if args.camera == "scientific":
        renderer.set_camera(ScientificCamera())
    else:
        renderer.set_camera(OrbitCamera())
    
    # Set transfer function based on shader
    if args.shader == "volume_grayscale":
        renderer.set_transfer_function(GrayscaleTransferFunction(args.field))
    else:
        renderer.set_transfer_function(ColorTransferFunction(args.field))
    
    # Render animation
    renderer.render_animation(
        args.output,
        fps=args.fps,
        duration=args.duration,
        speed=args.speed
    )


if __name__ == "__main__":
    main()
