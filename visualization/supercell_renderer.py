#!/usr/bin/env python3
"""
3D Volume Renderer for Supercell Simulations using ModernGL and volume ray marching.
"""
import moderngl as mgl
import moderngl_window as mglw
import numpy as np
import zarr
from pathlib import Path
import time
import math
from typing import Optional, Tuple

class SupercellVolumeRenderer(mglw.WindowConfig):
    """
    3D volume renderer for supercell thunderstorm simulations.
    Uses volume ray marching to render atmospheric fields.
    """
    title = "Supercell 3D Volume Renderer"
    gl_version = (3, 3)
    aspect_ratio = None
    vsync = True

    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument('--input', required=True, help='Zarr dataset path')
        parser.add_argument('--field', default='theta',
                          choices=['theta', 'qr', 'qv', 'qc', 'qi', 'qs', 'qh', 'u', 'v', 'w',
                                  'reflectivity_dbz', 'buoyancy', 'helicity'],
                          help='Primary field to visualize')
        parser.add_argument('--speed', type=float, default=1.0, help='Animation speed multiplier')
        parser.add_argument('--camera-distance', type=float, default=50.0, help='Initial camera distance')
        parser.add_argument('--list-fields', action='store_true', help='List available fields and exit')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Parse arguments
        self.zarr_path = Path(self.argv.input)
        self.primary_field = getattr(self.argv, 'field', 'theta')
        self.animation_speed = getattr(self.argv, 'speed', 1.0)
        self.camera_distance = getattr(self.argv, 'camera_distance', 50.0)
        self.list_fields_only = getattr(self.argv, 'list_fields', False)

        # Handle field listing request
        if self.list_fields_only:
            self.list_available_fields()
            # Exit after listing fields
            import sys
            sys.exit(0)

        # Load Zarr data
        self.load_zarr_data()

        # Rendering state
        self.current_frame = 0
        self.is_playing = True
        self.show_wireframe = False

        # Camera state
        self.camera_angle = 0.0
        self.camera_height = 20.0
        self.camera_target = np.array([0.0, 0.0, 10.0])

        # Transfer function parameters
        self.opacity_scale = 0.5  # Increased default opacity
        self.brightness = 2.0     # Increased default brightness

        # Initialize rendering
        self.setup_rendering()

        # Ensure all data is ready before starting
        print("Renderer initialization complete")
        print("=== CONTROLS ===")
        print("SPACE: Play/Pause | R: Reset | W: Wireframe")
        print("Up/Down: Opacity | Left/Right: Brightness | Mouse: Orbit")
        print("1-4: Presets | I: Help")
        print(f"Started with: Opacity={self.opacity_scale:.3f}, Brightness={self.brightness:.2f}")
        print("Try pressing '2' for high visibility preset!")

    def list_available_fields(self):
        """List all available fields in the Zarr dataset and exit"""
        try:
            store = zarr.open(str(self.zarr_path), mode='r')
            available_fields = []

            print(f"\nAvailable fields in {self.zarr_path}:")
            print("=" * 50)

            # Categorize fields
            thermo_fields = []
            moisture_fields = []
            wind_fields = []
            diagnostic_fields = []

            for field_name in store.keys():
                try:
                    field_data = store[field_name]
                    if hasattr(field_data, 'shape'):
                        shape = field_data.shape
                        dims = f"{shape}"
                    else:
                        # Group structure
                        keys = list(field_data.keys())
                        if keys:
                            shape = field_data[keys[0]].shape
                            dims = f"{len(keys)}×{shape}"
                        else:
                            dims = "empty"

                    # Categorize fields
                    if field_name in ['theta', 'theta_e', 'temperature']:
                        thermo_fields.append(f"  {field_name:<20} {dims}")
                    elif field_name in ['qv', 'qc', 'qr', 'qi', 'qs', 'qh']:
                        moisture_fields.append(f"  {field_name:<20} {dims}")
                    elif field_name in ['u', 'v', 'w']:
                        wind_fields.append(f"  {field_name:<20} {dims}")
                    elif field_name in ['reflectivity_dbz', 'buoyancy', 'helicity', 'srh', 'cape', 'cin']:
                        diagnostic_fields.append(f"  {field_name:<20} {dims}")
                    else:
                        available_fields.append(f"  {field_name:<20} {dims}")

                except Exception as e:
                    available_fields.append(f"  {field_name:<20} (error: {e})")

            if thermo_fields:
                print("Thermodynamic fields:")
                for field in thermo_fields:
                    print(field)

            if moisture_fields:
                print("\nMoisture/Precipitation fields:")
                for field in moisture_fields:
                    print(field)

            if wind_fields:
                print("\nWind fields:")
                for field in wind_fields:
                    print(field)

            if diagnostic_fields:
                print("\nDiagnostic fields:")
                for field in diagnostic_fields:
                    print(field)

            if available_fields:
                print("\nOther fields:")
                for field in available_fields:
                    print(field)

            print("\n" + "=" * 50)
            print("Use --field FIELD_NAME to visualize a specific field")
            print("Examples:")
            print("  --field theta      # Potential temperature")
            print("  --field qr         # Rain water mixing ratio")
            print("  --field qv         # Water vapor mixing ratio")

        except Exception as e:
            print(f"Error reading Zarr file: {e}")
            return

    def load_zarr_data(self):
        """Load simulation data from Zarr store"""
        if not self.zarr_path.exists():
            raise FileNotFoundError(f"Zarr file not found: {self.zarr_path}")

        self.store = zarr.open(str(self.zarr_path), mode='r')

        # Check if data exists
        if self.primary_field not in self.store:
            available_fields = list(self.store.keys())
            raise ValueError(f"Field '{self.primary_field}' not found. Available: {available_fields}")

        # Load dimensions - handle different zarr structures
        field_data = self.store[self.primary_field]

        if hasattr(field_data, 'shape'):
            # New structure: single 4D array (time, z, y, x)
            self.nt, self.nz, self.ny, self.nx = field_data.shape
            print(f"Loaded data: {self.nt} timesteps, {self.nx}x{self.ny}x{self.nz} grid (4D array)")
        elif hasattr(field_data, 'keys'):
            # Old structure: group with individual timestep arrays
            timestep_keys = sorted(field_data.keys(), key=lambda x: int(x))
            if not timestep_keys:
                raise ValueError(f"No timesteps found in {self.primary_field}")

            # Load first timestep to get dimensions
            first_timestep = field_data[timestep_keys[0]]
            self.nt = len(timestep_keys)
            _, self.nz, self.ny, self.nx = first_timestep.shape
            print(f"Loaded data: {self.nt} timesteps, {self.nx}x{self.ny}x{self.nz} grid (group structure)")
        else:
            raise ValueError(f"Unsupported zarr structure for {self.primary_field}")

        # Load coordinate arrays if available
        if 'x' in self.store:
            self.x_coords = self.store['x'][:]
            self.y_coords = self.store['y'][:]
            z_coords_full = self.store['z'][:]

            # Handle coordinate dimension mismatch
            if len(z_coords_full) != self.nz:
                print(f"Warning: z coordinate dimension mismatch. Data has {self.nz} levels, coords have {len(z_coords_full)} values.")
                # Truncate to match data dimensions (simpler approach)
                if len(z_coords_full) > self.nz:
                    # Truncate to match data dimensions
                    self.z_coords = z_coords_full[:self.nz]
                else:
                    # Pad with zeros if coords are too short (rare case)
                    padding = np.zeros(self.nz - len(z_coords_full))
                    self.z_coords = np.concatenate([z_coords_full, padding])
                print(f"Adjusted z coordinates to {len(self.z_coords)} values")
            else:
                self.z_coords = z_coords_full
        else:
            # Create default coordinates
            self.x_coords = np.linspace(-50, 50, self.nx)
            self.y_coords = np.linspace(-50, 50, self.ny)
            self.z_coords = np.linspace(0, 15, self.nz)

        # Pre-load some data for performance
        self.cache_data()

        # Find valid timesteps (not all NaN)
        self.find_valid_timesteps()
        print(f"Found {len(self.valid_timesteps)} valid timesteps out of {self.nt} total")

        # Pre-compute all transfer functions for smooth animation
        self.precompute_transfer_functions()

    def cache_data(self):
        """Cache frequently accessed data"""
        # Load primary field - handle different structures
        primary_field_data = self.store[self.primary_field]
        if hasattr(primary_field_data, 'shape'):
            # New structure: direct array access
            self.primary_data = primary_field_data[:]
        else:
            # Old structure: load from group
            timestep_keys = sorted(primary_field_data.keys(), key=lambda x: int(x))
            timestep_data = []
            for key in timestep_keys:
                data = primary_field_data[key][:]
                timestep_data.append(data)
            self.primary_data = np.stack(timestep_data, axis=0)

        # Load secondary fields for transfer function
        # Available fields from various physics modules
        available_fields = ['qr', 'qv', 'qc', 'qh', 'qi', 'qs', 'u', 'v', 'w']

        # Add diagnostic fields if available
        diagnostic_fields = ['reflectivity_dbz', 'buoyancy', 'helicity']
        available_fields.extend(diagnostic_fields)

        self.secondary_fields = {}
        for field_name in available_fields:
            if field_name in self.store:
                try:
                    field_data = self.store[field_name]
                    if hasattr(field_data, 'shape'):
                        # New structure: direct array access
                        self.secondary_fields[field_name] = field_data[:]
                        print(f"Loaded field: {field_name} {field_data.shape}")
                    else:
                        # Old structure: load from group
                        timestep_keys = sorted(field_data.keys(), key=lambda x: int(x))
                        timestep_data = []
                        for key in timestep_keys:
                            data = field_data[key][:]
                            timestep_data.append(data)
                        if timestep_data:
                            self.secondary_fields[field_name] = np.stack(timestep_data, axis=0)
                            print(f"Loaded field: {field_name} {len(timestep_data)} timesteps")
                except Exception as e:
                    print(f"Warning: Could not load field {field_name}: {e}")
                    continue

        # Report loaded fields
        if self.secondary_fields:
            loaded_fields = list(self.secondary_fields.keys())
            print(f"Loaded secondary fields: {', '.join(loaded_fields)}")
        else:
            print("No secondary fields found - visualization will use primary field only")

    def find_valid_timesteps(self):
        """Find timesteps that have valid (non-NaN) data"""
        self.valid_timesteps = []
        for t in range(self.nt):
            primary = self.primary_data[t]
            nan_fraction = np.isnan(primary).mean()
            # Consider a timestep valid if less than 50% of values are NaN
            if nan_fraction < 0.5:
                self.valid_timesteps.append(t)
        # If no valid timesteps found, use first timestep
        if not self.valid_timesteps:
            self.valid_timesteps = [0]

    def precompute_transfer_functions(self):
        """Pre-compute transfer functions for all valid timesteps to ensure smooth animation"""
        print(f"Pre-computing transfer functions for {len(self.valid_timesteps)} timesteps...")
        self.precomputed_rgba = {}

        for i, timestep in enumerate(self.valid_timesteps):
            if i % 10 == 0:  # Progress indicator
                print(f"  Processing timestep {i+1}/{len(self.valid_timesteps)} (t={timestep})")

            rgba_data = self.create_transfer_function(timestep)
            self.precomputed_rgba[timestep] = rgba_data

        print("Transfer function pre-computation complete")

    def setup_rendering(self):
        """Initialize OpenGL rendering pipeline"""
        # Create 3D volume texture
        self.volume_texture = self.ctx.texture3d(
            (self.nx, self.ny, self.nz), 4, dtype='f4'
        )
        self.volume_texture.filter = (mgl.LINEAR, mgl.LINEAR)
        self.volume_texture.repeat_x = False
        self.volume_texture.repeat_y = False
        self.volume_texture.repeat_z = False

        # Load shaders
        shader_dir = Path(__file__).parent / "shaders"
        vertex_shader = self.load_shader_file(shader_dir / "volume.vert")
        fragment_shader = self.load_shader_file(shader_dir / "volume.frag")

        # Create shader program
        self.volume_prog = self.ctx.program(
            vertex_shader=vertex_shader,
            fragment_shader=fragment_shader
        )

        # Create cube VAO for volume rendering
        self.cube_vao = self.create_cube_vao()

        # Create wireframe cube for debugging
        self.wireframe_vao = self.create_wireframe_cube_vao()

        # Initialize volume data
        self.update_volume_texture(0)

    def load_shader_file(self, shader_path: Path) -> str:
        """Load shader source from file"""
        if shader_path.exists():
            return shader_path.read_text()
        else:
            # Fallback to embedded shaders
            return self.get_embedded_shader(shader_path.name)

    def get_embedded_shader(self, filename: str) -> str:
        """Fallback embedded shaders"""
        if filename == "volume.vert":
            return """
            #version 330
            in vec3 in_position;
            out vec3 worldPos;
            out vec3 viewDir;

            uniform mat4 mvpMatrix;
            uniform vec3 cameraPos;

            void main() {
                worldPos = in_position;
                viewDir = normalize(in_position - cameraPos);
                gl_Position = mvpMatrix * vec4(in_position, 1.0);
            }
            """
        elif filename == "volume.frag":
            return """
            #version 330
            in vec3 worldPos;
            in vec3 viewDir;
            out vec4 fragColor;

            uniform sampler3D volumeTexture;
            uniform vec3 cameraPos;
            uniform float opacityScale;
            uniform float brightness;

            void main() {
                vec3 rayDir = normalize(viewDir);
                vec3 pos = cameraPos;

                vec4 accumulatedColor = vec4(0.0);
                float accumulatedAlpha = 0.0;

                // Ray march through volume
                for(int i = 0; i < 128; i++) {
                    // Transform to texture coordinates (volume is in [-1,1] cube)
                    vec3 texCoord = (pos + 1.0) * 0.5;

                    if(all(greaterThanEqual(texCoord, vec3(0.0))) &&
                       all(lessThanEqual(texCoord, vec3(1.0)))) {

                        vec4 sampleColor = texture(volumeTexture, texCoord);

                        // Apply transfer function
                        sampleColor.rgb *= brightness;
                        sampleColor.a *= opacityScale;

                        // Front-to-back blending
                        accumulatedColor += sampleColor * sampleColor.a * (1.0 - accumulatedAlpha);
                        accumulatedAlpha += sampleColor.a;

                        if(accumulatedAlpha > 0.95) break;
                    }

                    pos += rayDir * 0.02;  // Step size
                }

                fragColor = accumulatedColor;
            }
            """
        else:
            raise ValueError(f"Unknown shader: {filename}")

    def create_cube_vao(self):
        """Create VAO for rendering the volume cube"""
        # Cube vertices (front faces for proper ray marching)
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

    def create_wireframe_cube_vao(self):
        """Create wireframe cube for debugging"""
        # Wireframe edges
        vertices = np.array([
            -1, -1, -1,   1, -1, -1,
             1, -1, -1,   1,  1, -1,
             1,  1, -1,  -1,  1, -1,
            -1,  1, -1,  -1, -1, -1,
            -1, -1,  1,   1, -1,  1,
             1, -1,  1,   1,  1,  1,
             1,  1,  1,  -1,  1,  1,
            -1,  1,  1,  -1, -1,  1,
            -1, -1, -1,  -1, -1,  1,
             1, -1, -1,   1, -1,  1,
             1,  1, -1,   1,  1,  1,
            -1,  1, -1,  -1,  1,  1,
        ], dtype=np.float32)

        vbo = self.ctx.buffer(vertices.tobytes())
        return self.ctx.vertex_array(
            self.ctx.program(
                vertex_shader="""
                #version 330
                in vec3 in_position;
                uniform mat4 mvpMatrix;
                void main() {
                    gl_Position = mvpMatrix * vec4(in_position, 1.0);
                }
                """,
                fragment_shader="""
                #version 330
                out vec4 fragColor;
                void main() {
                    fragColor = vec4(1.0, 1.0, 1.0, 1.0);
                }
                """
            ),
            [(vbo, '3f', 'in_position')]
        )

    def create_transfer_function(self, timestep: int):
        """Create RGBA volume texture from physical fields"""
        # Get primary field data
        primary = self.primary_data[timestep]

        # Handle NaN values first
        nan_mask = np.isnan(primary)
        if np.any(nan_mask):
            print(f"Warning: timestep {timestep} has {np.sum(nan_mask)} NaN values")
            # Replace NaN values with a reasonable default
            primary = np.where(nan_mask, 300.0, primary)  # Use 300K as default temperature

        # Normalize primary field
        primary_min, primary_max = np.nanmin(primary), np.nanmax(primary)

        # Check for data corruption and apply field-specific normalization
        if self.primary_field == 'theta':
            # Theta should be reasonable temperature range (potential temperature)
            if primary_max > 1000 or primary_min < 200:
                print(f"Warning: timestep {timestep} has suspicious theta range: {primary_min:.1f} to {primary_max:.1f}")
                # Clamp to reasonable temperature range for potential temperature
                primary = np.clip(primary, 250, 500)
                # Re-compute min/max after clamping
                primary_min, primary_max = np.nanmin(primary), np.nanmax(primary)

            # Standard linear normalization for temperature
            if primary_max > primary_min:
                primary_norm = (primary - primary_min) / (primary_max - primary_min)
            else:
                primary_norm = np.zeros_like(primary)

        elif self.primary_field == 'reflectivity_dbz':
            # Radar reflectivity should be reasonable dBZ range (-30 to 60+)
            if primary_max > 80 or primary_min < -50:
                print(f"Warning: timestep {timestep} has suspicious reflectivity range: {primary_min:.1f} to {primary_max:.1f} dBZ")
                # Clamp to reasonable reflectivity range
                primary = np.clip(primary, -30, 60)
                # Re-compute min/max after clamping
                primary_min, primary_max = np.nanmin(primary), np.nanmax(primary)

            # For radar reflectivity, use a non-linear mapping to emphasize precipitation
            # Map -30 dBZ (noise floor) to 0, and higher reflectivity to 1
            # Use logarithmic-like scaling to show dynamic range
            reflectivity_threshold = -10.0  # dBZ, below this is mostly noise
            primary_norm = np.clip((primary - reflectivity_threshold) / (primary_max - reflectivity_threshold), 0.0, 1.0)

        else:
            # Default linear normalization for other fields
            if primary_max > primary_min:
                primary_norm = (primary - primary_min) / (primary_max - primary_min)
            else:
                primary_norm = np.zeros_like(primary)

        # Initialize RGBA channels
        r = primary_norm.copy()  # Red channel for primary field
        g = np.zeros_like(primary_norm)  # Green for secondary fields
        b = np.zeros_like(primary_norm)  # Blue for wind/turbulence
        a = np.zeros_like(primary_norm)  # Alpha for opacity

        # Debug: Check RGBA stats (only for first timestep)
        if timestep == self.valid_timesteps[0]:
            print(f"RGBA stats: R:{r.mean():.3f}, G:{g.mean():.3f}, B:{b.mean():.3f}, A:{a.mean():.3f}")

        # Add condensate fields (clouds, rain, ice, etc.)
        condensate_opacity = np.zeros_like(a)

        # Rain water mixing ratio
        if 'qr' in self.secondary_fields:
            qr = self.secondary_fields['qr'][timestep]
            qr_norm = np.clip(qr * 1000.0, 0.0, 1.0)  # Scale for visibility
            g += qr_norm * 0.6  # Green for precipitation
            condensate_opacity += qr_norm * 0.5

        # Cloud water mixing ratio
        if 'qc' in self.secondary_fields:
            qc = self.secondary_fields['qc'][timestep]
            qc_norm = np.clip(qc * 1000.0, 0.0, 1.0)
            g += qc_norm * 0.4  # Add to green channel for clouds
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

        # Water vapor field (blue channel secondary)
        if 'qv' in self.secondary_fields:
            qv = self.secondary_fields['qv'][timestep]
            # Water vapor is typically 5-20 g/kg, scale accordingly
            qv_norm = np.clip((qv - 0.005) / 0.015, 0.0, 1.0)  # Normalize to reasonable range
            b += qv_norm * 0.5  # Add to blue channel

        # Wind field visualization (primary blue channel)
        if 'u' in self.secondary_fields and 'v' in self.secondary_fields and 'w' in self.secondary_fields:
            u = self.secondary_fields['u'][timestep]
            v = self.secondary_fields['v'][timestep]
            w = self.secondary_fields['w'][timestep]

            wind_speed = np.sqrt(u**2 + v**2 + w**2)
            wind_norm = np.clip(wind_speed / 50.0, 0.0, 1.0)  # Scale for typical storm winds

            b = np.maximum(b, wind_norm)  # Blue for wind intensity

            # Add some wind contribution to opacity
            condensate_opacity = np.maximum(condensate_opacity, wind_norm * 0.2)

        # Add condensate opacity to final alpha
        a = np.maximum(a, condensate_opacity)

        # Ensure minimum opacity for structure visibility
        a = np.maximum(a, 0.1)  # Increased minimum opacity

        # Stack into RGBA
        rgba = np.stack([r, g, b, a], axis=-1)

        return rgba.astype(np.float32)

    def update_volume_texture(self, timestep: int):
        """Update 3D texture with current timestep data"""
        # Use precomputed data for smooth animation
        if timestep in self.precomputed_rgba:
            rgba_data = self.precomputed_rgba[timestep]
        else:
            # Fallback to on-demand computation if not precomputed
            rgba_data = self.create_transfer_function(timestep)
        self.volume_texture.write(rgba_data.tobytes())

    def create_view_matrix(self) -> np.ndarray:
        """Create view matrix for orbiting camera"""
        # Camera position (orbiting around storm)
        eye_x = self.camera_distance * math.cos(self.camera_angle)
        eye_y = self.camera_distance * math.sin(self.camera_angle)
        eye_z = self.camera_height

        eye = np.array([eye_x, eye_y, eye_z])
        center = self.camera_target
        up = np.array([0.0, 0.0, 1.0])

        return self.look_at(eye, center, up)

    def create_projection_matrix(self) -> np.ndarray:
        """Create perspective projection matrix"""
        fov = math.radians(45.0)
        aspect = self.wnd.size[0] / self.wnd.size[1]
        near = 0.1
        far = 200.0

        return self.perspective(fov, aspect, near, far)

    @staticmethod
    def look_at(eye: np.ndarray, center: np.ndarray, up: np.ndarray) -> np.ndarray:
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

    @staticmethod
    def perspective(fov: float, aspect: float, near: float, far: float) -> np.ndarray:
        """Create perspective projection matrix"""
        f = 1.0 / math.tan(fov / 2.0)
        m = np.zeros((4, 4))

        m[0, 0] = f / aspect
        m[1, 1] = f
        m[2, 2] = (far + near) / (near - far)
        m[2, 3] = (2.0 * far * near) / (near - far)
        m[3, 2] = -1.0

        return m.astype(np.float32)

    def on_render(self, time: float, frame_time: float):
        """Main render function"""
        self.ctx.clear(0.1, 0.1, 0.2, 1.0)  # Dark blue background
        self.ctx.enable(mgl.DEPTH_TEST)
        self.ctx.enable(mgl.BLEND)
        self.ctx.blend_func = mgl.SRC_ALPHA, mgl.ONE_MINUS_SRC_ALPHA

        # Update animation
        if self.is_playing and self.valid_timesteps:
            frame_index = int(time * 30 * self.animation_speed) % len(self.valid_timesteps)
            self.current_frame = self.valid_timesteps[frame_index]
            self.update_volume_texture(self.current_frame)

            # Slowly orbit camera
            self.camera_angle += frame_time * 0.2

        # Create transformation matrices
        view_matrix = self.create_view_matrix()
        proj_matrix = self.create_projection_matrix()
        mvp_matrix = proj_matrix @ view_matrix

        # Set shader uniforms
        self.volume_prog['mvpMatrix'].write(mvp_matrix.tobytes())
        self.volume_prog['cameraPos'].write(view_matrix[:3, 3].tobytes())
        self.volume_prog['opacityScale'].value = self.opacity_scale
        self.volume_prog['brightness'].value = self.brightness
        self.volume_texture.use(0)

        # Render volume
        self.cube_vao.render(mgl.TRIANGLES)

        # Render wireframe for debugging
        if self.show_wireframe:
            self.wireframe_vao.render(mgl.LINES)

    def on_key_event(self, key, action, modifiers):
        """Handle keyboard input"""
        if action == self.wnd.keys.ACTION_PRESS:
            if key == self.wnd.keys.SPACE:
                self.is_playing = not self.is_playing
                print(f"Animation: {'Playing' if self.is_playing else 'Paused'}")
            elif key == self.wnd.keys.R:
                self.current_frame = 0
                self.camera_angle = 0.0
                print("Reset camera and animation")
            elif key == self.wnd.keys.W:
                self.show_wireframe = not self.show_wireframe
                print(f"Wireframe: {'ON' if self.show_wireframe else 'OFF'}")
            elif key == self.wnd.keys.UP:
                self.opacity_scale *= 1.2
                print(f"Opacity: {self.opacity_scale:.3f}")
            elif key == self.wnd.keys.DOWN:
                self.opacity_scale /= 1.2
                print(f"Opacity: {self.opacity_scale:.3f}")
            elif key == self.wnd.keys.LEFT:
                self.brightness *= 1.1
                print(f"Brightness: {self.brightness:.2f}")
            elif key == self.wnd.keys.RIGHT:
                self.brightness /= 1.1
                print(f"Brightness: {self.brightness:.2f}")
            # Number keys for presets
            elif key == self.wnd.keys._1:
                self.opacity_scale = 0.1
                self.brightness = 1.0
                print("Preset 1: Default settings")
            elif key == self.wnd.keys._2:
                self.opacity_scale = 0.5
                self.brightness = 2.0
                print("Preset 2: High visibility")
            elif key == self.wnd.keys._3:
                self.opacity_scale = 1.0
                self.brightness = 3.0
                print("Preset 3: Maximum visibility")
            elif key == self.wnd.keys._4:
                self.opacity_scale = 0.01
                self.brightness = 0.5
                print("Preset 4: Low visibility (debug)")
            elif key == self.wnd.keys.I:
                print("=== CONTROLS ===")
                print("SPACE: Play/Pause animation")
                print("R: Reset camera/animation")
                print("W: Toggle wireframe")
                print("↑↓: Adjust opacity")
                print("Left/Right arrows: Adjust brightness")
                print("1-4: Preset settings")
                print("I: Show this help")
                print("Mouse: Orbit camera")
                print(f"Current: Opacity={self.opacity_scale:.3f}, Brightness={self.brightness:.2f}")

    def on_mouse_drag_event(self, x, y, dx, dy):
        """Handle mouse drag for camera control"""
        self.camera_angle += dx * 0.01
        self.camera_height += dy * 0.1
        self.camera_height = np.clip(self.camera_height, 5.0, 100.0)

def main():
    """Main entry point"""
    mglw.run_window_config(SupercellVolumeRenderer)

if __name__ == "__main__":
    main()
