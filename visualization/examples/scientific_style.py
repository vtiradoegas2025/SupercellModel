"""
Scientific-style visualization example.

Creates visualizations similar to Lewellen et al. (2008) with:
- Grayscale volume rendering for downdraft column
- Dense contour plots for turbulence structures
- Side-view camera angles
"""

import argparse
import numpy as np
import moderngl as mgl
from pathlib import Path
from PIL import Image
import tempfile
import subprocess
import shutil
from tqdm import tqdm

from ..core.data_manager import DataManager
from ..core.render_engine import RenderEngine
from ..core.transfer_function import GrayscaleTransferFunction
from ..core.camera import ScientificCamera
from ..passes.volume_pass import VolumeRenderPass
from ..passes.contour_pass import ContourRenderPass


class ScientificStyleRenderer:
    """
    Renderer for scientific-style visualizations.
    """
    
    def __init__(self, data_path: str, width: int = 1920, height: int = 1080):
        """
        Initialize scientific style renderer.
        
        Args:
            data_path: Path to simulation data
            width: Viewport width
            height: Viewport height
        """
        self.data_path = Path(data_path)
        self.width = width
        self.height = height
        
        # Initialize data manager
        self.data_manager = DataManager(data_path)
        
        # Initialize render engine
        shader_dir = Path(__file__).parent.parent / "shaders"
        self.engine = RenderEngine(width, height, shader_dir=shader_dir)
        
        # Setup camera (side view for downdraft)
        self.camera = ScientificCamera()
        self.camera.set_side_view()
        self.engine.set_camera(self.camera)
        
        # Setup volume texture
        nx, ny, nz = self.data_manager.get_grid_dims()
        self.volume_texture = self.engine.ctx.texture3d(
            (nx, ny, nz), 4, dtype='f4'
        )
        self.volume_texture.filter = (mgl.LINEAR, mgl.LINEAR)
        
        # Setup transfer function (grayscale)
        self.transfer_function = GrayscaleTransferFunction("theta")
        
        # Setup volume render pass (grayscale)
        self.volume_pass = VolumeRenderPass(
            self.engine.shader_manager,
            self.transfer_function,
            self.volume_texture,
            vertex_shader="volume",
            fragment_shader="volume_grayscale",
            opacity_scale=0.3,
            brightness=1.5
        )
        self.engine.add_render_pass(self.volume_pass)
        
        # Setup contour pass for turbulence
        # Use a separate texture for contour data
        self.contour_texture = self.engine.ctx.texture3d(
            (nx, ny, nz), 4, dtype='f4'
        )
        self.contour_texture.filter = (mgl.LINEAR, mgl.LINEAR)
        
        self.contour_pass = ContourRenderPass(
            self.engine.shader_manager,
            self.contour_texture,
            contour_spacing=0.05,
            contour_width=0.01,
            contour_color=(0.0, 0.0, 0.0),  # Black contours
            opacity=0.8
        )
        # Note: Contour pass would need to be rendered separately or composited
        # For now, we'll render volume only
    
    def render_frame(self, timestep: int = 0, field_name: str = "theta") -> np.ndarray:
        """
        Render a single frame.
        
        Args:
            timestep: Timestep index
            field_name: Field to visualize
            
        Returns:
            Rendered image as numpy array
        """
        # Load field data
        field_data = self.data_manager.get_field(field_name, timestep)
        
        # Apply transfer function
        rgba_data = self.transfer_function.apply(field_data)
        
        # Update volume texture
        self.volume_texture.write(rgba_data.tobytes())
        
        # Render
        pixels = self.engine.render_frame(timestep=timestep)
        
        return pixels
    
    def render_animation(self, output_path: str, fps: int = 10, max_timesteps: int = 30):
        """
        Render animation with scientific style.
        
        Args:
            output_path: Output video file path
            fps: Frame rate
            max_timesteps: Maximum number of timesteps to render
        """
        # Determine number of timesteps
        if self.data_manager.format == 'zarr':
            sample_field = self.data_manager.get_all_timesteps("theta")
            num_timesteps = min(sample_field.shape[0], max_timesteps) if sample_field is not None else 1
        else:
            import glob
            step_dirs = glob.glob(str(self.data_path / "step_*"))
            num_timesteps = min(len(step_dirs), max_timesteps)
        
        # Create temporary directory for frames
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            print(f"Rendering {num_timesteps} frames with scientific style...")
            
            # Render each frame
            for timestep in tqdm(range(num_timesteps), desc="Rendering"):
                pixels = self.render_frame(timestep)
                
                # Save frame
                img = Image.fromarray(pixels, 'RGBA')
                frame_path = temp_path / f"frame_{timestep:06d}.png"
                img.save(frame_path)
            
            print("Encoding video...")
            
            # Encode with FFmpeg
            ffmpeg_path = shutil.which("ffmpeg")
            if ffmpeg_path is None:
                raise RuntimeError("FFmpeg not found. Install FFmpeg to create videos.")
            
            cmd = [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(temp_path / "frame_%06d.png"),
                "-c:v", "libx264",
                "-preset", "slow",
                "-crf", "22",
                "-pix_fmt", "yuv420p",
                str(output_path)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"FFmpeg error: {result.stderr}")
                raise RuntimeError("Video encoding failed")
        
        print(f"Scientific-style animation saved to: {output_path}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(description="Create scientific-style visualizations")
    parser.add_argument("--input", required=True, help="Input data path (Zarr or NPY)")
    parser.add_argument("--output", required=True, help="Output video file")
    parser.add_argument("--fps", type=int, default=10, help="Frame rate")
    parser.add_argument("--max-timesteps", type=int, default=30, help="Max timesteps")
    parser.add_argument("--width", type=int, default=1920, help="Width")
    parser.add_argument("--height", type=int, default=1080, help="Height")
    
    args = parser.parse_args()
    
    renderer = ScientificStyleRenderer(
        args.input,
        width=args.width,
        height=args.height
    )
    
    renderer.render_animation(
        args.output,
        fps=args.fps,
        max_timesteps=args.max_timesteps
    )


if __name__ == "__main__":
    main()
