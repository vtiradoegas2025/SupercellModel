#!/usr/bin/env python3
"""
Offline 3D volume renderer for supercell simulations.
Renders high-quality videos from Zarr datasets using OpenGL volume ray marching.

This script uses the new modular engine internally while maintaining backward compatibility.
All existing CLI interface and behavior remain unchanged.
"""
import argparse
import os
import subprocess
import shutil
import tempfile
from pathlib import Path
import numpy as np
from PIL import Image

# Conditional imports for OpenGL dependencies
try:
    import moderngl as mgl
    import zarr
    from tqdm import tqdm
    OPENGL_AVAILABLE = True
except ImportError as e:
    OPENGL_AVAILABLE = False
    missing_module = str(e).split("'")[1] if "'" in str(e) else "unknown"
    # Create dummy classes
    class mgl:
        pass
    class zarr:
        pass
    class tqdm:
        def __init__(self, *args, **kwargs):
            pass
        def __enter__(self):
            return self
        def __exit__(self, *args):
            pass

# Import new engine components
ENGINE_AVAILABLE = False
try:
    # Try relative import first (when used as module)
    try:
        from .renderers.offline_renderer import OfflineRenderer as NewOfflineRenderer
        ENGINE_AVAILABLE = True
    except (ImportError, ModuleNotFoundError):
        # Fallback to absolute import (when run as script)
        try:
            import sys
            from pathlib import Path
            vis_path = Path(__file__).parent
            if str(vis_path) not in sys.path:
                sys.path.insert(0, str(vis_path))
            from renderers.offline_renderer import OfflineRenderer as NewOfflineRenderer
            ENGINE_AVAILABLE = True
        except (ImportError, ModuleNotFoundError):
            # Engine components not available
            ENGINE_AVAILABLE = False
except Exception:
    ENGINE_AVAILABLE = False

class OfflineVolumeRenderer:
    """
    Renders 3D volume animations to video files without display window.
    
    Uses new modular engine internally while maintaining same interface.
    """
    def __init__(self, zarr_path: str, width: int = 1920, height: int = 1080):
        self.zarr_path = Path(zarr_path)
        self.width = width
        self.height = height

        # Use new OfflineRenderer if available
        if ENGINE_AVAILABLE:
            try:
                self._new_renderer = NewOfflineRenderer(
                    str(zarr_path),
                    primary_field="theta",
                    shader="volume_color",
                    width=width,
                    height=height
                )
                # Map old interface to new renderer
                self.ctx = self._new_renderer.engine.ctx
                self.fbo = self._new_renderer.engine.fbo
                self.nt = self._new_renderer.num_timesteps
                print("Using new modular engine for rendering")
                return
            except Exception as e:
                print(f"Warning: New engine failed, using fallback: {e}")
                self._new_renderer = None
        
        # Fallback to original implementation
        self._new_renderer = None
        
        # Create headless OpenGL context
        self.ctx = mgl.create_standalone_context()
        self.ctx.enable(mgl.DEPTH_TEST)
        self.ctx.enable(mgl.BLEND)
        self.ctx.blend_func = mgl.SRC_ALPHA, mgl.ONE_MINUS_SRC_ALPHA

        # Load data
        self.load_data()

        # Create FBO for offscreen rendering
        self.fbo = self.ctx.framebuffer(
            color_attachments=[self.ctx.texture((width, height), 4)],
            depth_attachment=self.ctx.depth_texture((width, height))
        )

        # Setup rendering (similar to interactive renderer)
        self.setup_rendering()

    def load_data(self):
        """Load Zarr data"""
        if not self.zarr_path.exists():
            raise FileNotFoundError(f"Zarr file not found: {self.zarr_path}")
        
        self.store = zarr.open(str(self.zarr_path), mode='r')
        
        # Validate that required field exists
        if 'theta' not in self.store:
            available_fields = list(self.store.keys())
            raise ValueError(f"Required field 'theta' not found in Zarr file. Available fields: {available_fields}")
        
        # Validate data dimensions
        theta_shape = self.store['theta'].shape
        if len(theta_shape) != 4:
            raise ValueError(f"Expected 4D array (time, z, y, x), got shape {theta_shape}")
        
        self.nt, self.nz, self.ny, self.nx = theta_shape
        
        if self.nt == 0:
            raise ValueError("No timesteps found in data")
        if self.nx == 0 or self.ny == 0 or self.nz == 0:
            raise ValueError(f"Invalid grid dimensions: {self.nx}x{self.ny}x{self.nz}")
        
        self.primary_data = self.store['theta'][:]

        # Load secondary fields - expanded to include all atmospheric fields
        available_fields = ['qr', 'qv', 'qc', 'qh', 'qi', 'qs', 'u', 'v', 'w']
        diagnostic_fields = ['reflectivity_dbz', 'buoyancy', 'helicity']
        available_fields.extend(diagnostic_fields)

        self.secondary_fields = {}
        for field_name in available_fields:
            if field_name in self.store:
                try:
                    self.secondary_fields[field_name] = self.store[field_name][:]
                    print(f"Loaded field: {field_name}")
                except Exception as e:
                    print(f"Warning: Could not load field {field_name}: {e}")
                    continue

    def setup_rendering(self):
        """Setup rendering pipeline"""
        # Create volume texture
        self.volume_texture = self.ctx.texture3d(
            (self.nx, self.ny, self.nz), 4, dtype='f4'
        )
        self.volume_texture.filter = (mgl.LINEAR, mgl.LINEAR)

        # Load shaders - try new location first, then fallback
        shader_dir = Path(__file__).parent / "shaders"
        vert_path = shader_dir / "volume" / "volume.vert"
        frag_path = shader_dir / "volume" / "volume_color.frag"
        if not vert_path.exists():
            vert_path = shader_dir / "volume.vert"
        if not frag_path.exists():
            frag_path = shader_dir / "volume.frag"
        
        vertex_shader = vert_path.read_text()
        fragment_shader = frag_path.read_text()

        # Create shader program
        self.volume_prog = self.ctx.program(
            vertex_shader=vertex_shader,
            fragment_shader=fragment_shader
        )

        # Create cube VAO
        self.cube_vao = self.create_cube_vao()

    def create_cube_vao(self):
        """Create VAO for volume cube"""
        vertices = np.array([
            # Front face
            -1, -1,  1,   1, -1,  1,   1,  1,  1,
            -1, -1,  1,   1,  1,  1,  -1,  1,  1,
            # Right face
             1, -1,  1,   1, -1, -1,   1,  1, -1,
             1, -1,  1,   1,  1, -1,   1,  1,  1,
            # Back face
             1, -1, -1,  -1, -1, -1,  -1,  1, -1,
             1, -1, -1,  -1,  1, -1,   1,  1, -1,
            # Left face
            -1, -1, -1,  -1, -1,  1,  -1,  1,  1,
            -1, -1, -1,  -1,  1,  1,  -1,  1, -1,
            # Top face
            -1,  1,  1,   1,  1,  1,   1,  1, -1,
            -1,  1,  1,   1,  1, -1,  -1,  1, -1,
            # Bottom face
            -1, -1, -1,   1, -1, -1,   1, -1,  1,
            -1, -1, -1,   1, -1,  1,  -1, -1,  1,
        ], dtype=np.float32)

        vbo = self.ctx.buffer(vertices.tobytes())
        return self.ctx.vertex_array(
            self.volume_prog,
            [(vbo, '3f', 'in_position')]
        )

    def create_transfer_function(self, timestep: int):
        """Create RGBA volume data for rendering"""
        # Get primary field data
        primary = self.primary_data[timestep]

        # Normalize
        primary_min, primary_max = np.nanmin(primary), np.nanmax(primary)
        if primary_max > primary_min:
            primary_norm = (primary - primary_min) / (primary_max - primary_min)
        else:
            primary_norm = np.zeros_like(primary)

        # Initialize RGBA
        r = primary_norm.copy()
        g = np.zeros_like(primary_norm)
        b = np.zeros_like(primary_norm)
        a = np.zeros_like(primary_norm)

        # Add condensate fields (clouds, rain, ice, etc.)
        condensate_opacity = np.zeros_like(a)

        # Rain water mixing ratio
        if 'qr' in self.secondary_fields:
            qr = self.secondary_fields['qr'][timestep]
            qr_norm = np.clip(qr * 1000.0, 0.0, 1.0)
            g += qr_norm * 0.6  # Green for precipitation
            condensate_opacity += qr_norm * 0.5

        # Cloud water mixing ratio
        if 'qc' in self.secondary_fields:
            qc = self.secondary_fields['qc'][timestep]
            qc_norm = np.clip(qc * 1000.0, 0.0, 1.0)
            g += qc_norm * 0.4  # Add to green for clouds
            condensate_opacity += qc_norm * 0.3

        # Cloud ice mixing ratio
        if 'qi' in self.secondary_fields:
            qi = self.secondary_fields['qi'][timestep]
            qi_norm = np.clip(qi * 1000.0, 0.0, 1.0)
            g += qi_norm * 0.3  # Lighter green for ice
            condensate_opacity += qi_norm * 0.4

        # Snow mixing ratio
        if 'qs' in self.secondary_fields:
            qs = self.secondary_fields['qs'][timestep]
            qs_norm = np.clip(qs * 1000.0, 0.0, 1.0)
            g += qs_norm * 0.2  # Even lighter for snow
            condensate_opacity += qs_norm * 0.3

        # Hail/graupel mixing ratio
        if 'qh' in self.secondary_fields:
            qh = self.secondary_fields['qh'][timestep]
            qh_norm = np.clip(qh * 1000.0, 0.0, 1.0)
            g += qh_norm * 0.5  # Distinct green for hail
            condensate_opacity += qh_norm * 0.6

        # Water vapor field
        if 'qv' in self.secondary_fields:
            qv = self.secondary_fields['qv'][timestep]
            qv_norm = np.clip((qv - 0.005) / 0.015, 0.0, 1.0)
            b += qv_norm * 0.5  # Add to blue channel

        # Wind field visualization
        if 'u' in self.secondary_fields and 'v' in self.secondary_fields and 'w' in self.secondary_fields:
            u = self.secondary_fields['u'][timestep]
            v = self.secondary_fields['v'][timestep]
            w = self.secondary_fields['w'][timestep]

            wind_speed = np.sqrt(u**2 + v**2 + w**2)
            wind_norm = np.clip(wind_speed / 50.0, 0.0, 1.0)

            b = np.maximum(b, wind_norm)  # Blue for wind intensity
            condensate_opacity = np.maximum(condensate_opacity, wind_norm * 0.2)

        # Add condensate opacity to final alpha
        a = np.maximum(a, condensate_opacity)

        # Minimum opacity
        a = np.maximum(a, 0.05)

        return np.stack([r, g, b, a], axis=-1).astype(np.float32)

    def update_volume_texture(self, timestep: int):
        """Update volume texture with current timestep"""
        rgba_data = self.create_transfer_function(timestep)
        self.volume_texture.write(rgba_data.tobytes())

    def create_view_matrix(self, angle: float, distance: float = 50.0, height: float = 20.0):
        """Create orbiting view matrix"""
        import math
        eye_x = distance * math.cos(angle)
        eye_y = distance * math.sin(angle)
        eye_z = height

        eye = np.array([eye_x, eye_y, eye_z])
        center = np.array([0.0, 0.0, 10.0])
        up = np.array([0.0, 0.0, 1.0])

        return self.look_at(eye, center, up)

    def create_projection_matrix(self):
        """Create perspective projection"""
        import math
        fov = math.radians(45.0)
        aspect = self.width / self.height
        near = 0.1
        far = 200.0

        f = 1.0 / math.tan(fov / 2.0)
        m = np.zeros((4, 4))

        m[0, 0] = f / aspect
        m[1, 1] = f
        m[2, 2] = (far + near) / (near - far)
        m[2, 3] = (2.0 * far * near) / (near - far)
        m[3, 2] = -1.0

        return m.astype(np.float32)

    @staticmethod
    def look_at(eye, center, up):
        """Create look-at matrix"""
        f = center - eye
        f = f / np.linalg.norm(f)

        s = np.cross(f, up)
        s = s / np.linalg.norm(s)

        u = np.cross(s, f)

        m = np.eye(4)
        m[0, :3] = s
        m[1, :3] = u
        m[2, :3] = -f
        m[:3, 3] = -np.dot(np.array([s, u, -f]).T, eye)

        return m.astype(np.float32)

    def render_frame(self, timestep: int, camera_angle: float = 0.0):
        """Render a single frame to the FBO"""
        # Update volume data
        self.update_volume_texture(timestep)

        # Bind FBO
        self.fbo.use()

        # Clear
        self.ctx.clear(0.1, 0.1, 0.2, 1.0)

        # Create matrices
        view_matrix = self.create_view_matrix(camera_angle)
        proj_matrix = self.create_projection_matrix()
        mvp_matrix = proj_matrix @ view_matrix

        # Set uniforms
        self.volume_prog['mvpMatrix'].write(mvp_matrix.tobytes())
        self.volume_prog['cameraPos'].write(view_matrix[:3, 3].tobytes())
        self.volume_prog['opacityScale'].value = 0.1
        self.volume_prog['brightness'].value = 1.0
        self.volume_prog['time'].value = timestep * 0.1

        # Bind texture
        self.volume_texture.use(0)

        # Render
        self.cube_vao.render(mgl.TRIANGLES)

        # Read back pixels
        pixels = self.fbo.read(components=4, dtype='f4')
        pixels = pixels.reshape((self.height, self.width, 4))

        # Convert to 8-bit RGBA
        pixels = np.clip(pixels * 255, 0, 255).astype(np.uint8)

        return pixels

    def render_animation(self, output_path: str, fps: int = 30, duration: float = None,
                        speed: float = 1.0):
        """Render full animation to video"""
        # Use new renderer if available
        if self._new_renderer is not None:
            try:
                self._new_renderer.render_animation(
                    output_path,
                    fps=fps,
                    duration=duration,
                    speed=speed,
                    camera_rotation=True
                )
                return
            except Exception as e:
                print(f"Warning: New renderer failed, using fallback: {e}")
        
        # Fallback to original implementation
        if duration is None:
            # Render entire simulation
            frame_count = self.nt
        else:
            frame_count = int(duration * fps)

        # Create temporary directory for frames
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            print(f"Rendering {frame_count} frames...")

            # Render each frame with progress bar
            for frame_idx in tqdm(range(frame_count), desc="Rendering frames", unit="frame"):
                timestep = int(frame_idx * speed) % self.nt
                camera_angle = frame_idx * 0.02  # Slow camera rotation

                # Render frame
                pixels = self.render_frame(timestep, camera_angle)

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

            # Encode to video with ffmpeg
            cmd = [
                "ffmpeg", "-y",
                "-framerate", str(fps),
                "-i", str(temp_path / "frame_%06d.png"),
                "-c:v", "libx264",
                "-preset", "slow",
                "-crf", "22",
                "-pix_fmt", "yuv420p",
                "-vf", "scale=1920:1080:force_original_aspect_ratio=decrease,pad=1920:1080:(ow-iw)/2:(oh-ih)/2",
                str(output_path)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"FFmpeg error: {result.stderr}")
                raise RuntimeError("Video encoding failed")

        print(f"Video saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Render 3D supercell animations")
    parser.add_argument("--input", required=True, help="Zarr dataset path")
    parser.add_argument("--output", required=True, help="Output video file")
    parser.add_argument("--fps", type=int, default=30, help="Frame rate")
    parser.add_argument("--duration", type=float, help="Duration in seconds")
    parser.add_argument("--speed", type=float, default=1.0, help="Animation speed multiplier")
    parser.add_argument("--width", type=int, default=1920, help="Video width")
    parser.add_argument("--height", type=int, default=1080, help="Video height")

    args = parser.parse_args()

    # Check dependencies after parsing (allows --help to work)
    if not OPENGL_AVAILABLE:
        print("ERROR: OpenGL dependencies not installed!")
        print("Required packages: moderngl, zarr, tqdm")
        print("Install with: pip install moderngl zarr tqdm")
        print("\nFor a simple 2D animation without OpenGL, use:")
        print("  python visualization/create_animation.py --input <data_path> --field theta --output animation.gif")
        return 1

    # Create renderer
    renderer = OfflineVolumeRenderer(args.input, args.width, args.height)

    # Render animation
    renderer.render_animation(
        args.output,
        fps=args.fps,
        duration=args.duration,
        speed=args.speed
    )

if __name__ == "__main__":
    main()
