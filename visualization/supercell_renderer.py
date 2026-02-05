#!/usr/bin/env python3
"""
3D Volume Renderer for Supercell Simulations using ModernGL and volume ray marching.

This renderer uses the new modular engine internally while maintaining backward compatibility.
All existing features, controls, and CLI interface remain unchanged.
"""
# Conditional imports for OpenGL dependencies
try:
    import moderngl as mgl
    import moderngl_window as mglw
    import zarr
    OPENGL_AVAILABLE = True
except ImportError as e:
    OPENGL_AVAILABLE = False
    missing_module = str(e).split("'")[1] if "'" in str(e) else "unknown"
    # Create dummy classes for type hints when imports fail
    class mgl:
        pass
    class mglw:
        class WindowConfig:
            pass
    class zarr:
        pass

import numpy as np
from pathlib import Path
import time
import math
from typing import Optional, Tuple

# Import new engine components for internal use
ENGINE_AVAILABLE = False
try:
    # Try relative import first (when used as module)
    try:
        from .core.data_manager import DataManager
        from .core.shader_manager import ShaderManager
        from .core.camera import OrbitCamera
        ENGINE_AVAILABLE = True
    except (ImportError, ModuleNotFoundError):
        # Fallback to absolute import (when run as script)
        try:
            import sys
            from pathlib import Path
            vis_path = Path(__file__).parent
            if str(vis_path) not in sys.path:
                sys.path.insert(0, str(vis_path))
            from core.data_manager import DataManager
            from core.shader_manager import ShaderManager
            from core.camera import OrbitCamera
            ENGINE_AVAILABLE = True
        except (ImportError, ModuleNotFoundError):
            # Engine components not available (e.g., missing dependencies)
            ENGINE_AVAILABLE = False
except Exception:
    # Any other error - fallback to original implementation
    ENGINE_AVAILABLE = False

# Define base class conditionally
if OPENGL_AVAILABLE:
    _BaseClass = mglw.WindowConfig
else:
    # Dummy base class for help/imports to work
    class _BaseClass:
        pass

class SupercellVolumeRenderer(_BaseClass):
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
        # Check dependencies - will be caught by main() if help is requested
        if not OPENGL_AVAILABLE:
            print("ERROR: OpenGL dependencies not installed!")
            print("Required packages: moderngl, moderngl-window, zarr")
            print("Install with: pip install moderngl moderngl-window zarr")
            print("\nFor a simple 2D animation without OpenGL, use:")
            print("  python visualization/create_animation.py --input <data_path> --field theta --output animation.gif")
            import sys
            sys.exit(1)
        
        super().__init__(**kwargs)

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
        
        # Initialize camera system if available
        if ENGINE_AVAILABLE:
            try:
                self._orbit_camera = OrbitCamera(
                    distance=self.camera_distance,
                    height=self.camera_height,
                    angle=self.camera_angle,
                    target=self.camera_target
                )
            except:
                self._orbit_camera = None
        else:
            self._orbit_camera = None

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
        # Use DataManager if available for source-of-truth alignment
        if ENGINE_AVAILABLE:
            try:
                data_mgr = DataManager(str(self.zarr_path))
                available_fields = data_mgr.list_available_fields()
                
                print(f"\nAvailable fields in {self.zarr_path}:")
                print("=" * 50)
                
                # Categorize fields (same as original)
                thermo_fields = []
                moisture_fields = []
                wind_fields = []
                diagnostic_fields = []
                other_fields = []
                
                for field_name in available_fields:
                    try:
                        field_data = data_mgr.get_all_timesteps(field_name)
                        if field_data is not None:
                            dims = f"{field_data.shape}"
                        else:
                            dims = "unknown"
                        
                        # Categorize
                        if field_name in ['theta', 'theta_e', 'temperature']:
                            thermo_fields.append(f"  {field_name:<20} {dims}")
                        elif field_name in ['qv', 'qc', 'qr', 'qi', 'qs', 'qh', 'qg']:
                            moisture_fields.append(f"  {field_name:<20} {dims}")
                        elif field_name in ['u', 'v', 'w']:
                            wind_fields.append(f"  {field_name:<20} {dims}")
                        elif field_name in ['reflectivity_dbz', 'radar', 'buoyancy', 'helicity', 'srh', 'cape', 'cin']:
                            diagnostic_fields.append(f"  {field_name:<20} {dims}")
                        else:
                            other_fields.append(f"  {field_name:<20} {dims}")
                    except Exception as e:
                        other_fields.append(f"  {field_name:<20} (error: {e})")
                
                # Print categorized
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
                if other_fields:
                    print("\nOther fields:")
                    for field in other_fields:
                        print(field)
                
                print("\n" + "=" * 50)
                print("Use --field FIELD_NAME to visualize a specific field")
                print("Examples:")
                print("  --field theta      # Potential temperature")
                print("  --field qr         # Rain water mixing ratio")
                print("  --field qv         # Water vapor mixing ratio")
                return
            except Exception as e:
                print(f"Warning: DataManager field listing failed: {e}, using fallback")
        
        # Fallback to original implementation
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
        # Use DataManager for source-of-truth alignment
        use_fallback = False
        if ENGINE_AVAILABLE:
            try:
                self.data_manager = DataManager(str(self.zarr_path))
                self.store = self.data_manager.store if hasattr(self.data_manager, 'store') else zarr.open(str(self.zarr_path), mode='r')
                
                # Get dimensions from DataManager
                self.nx, self.ny, self.nz = self.data_manager.get_grid_dims()
                
                # Get number of timesteps
                all_timesteps = self.data_manager.get_all_timesteps(self.primary_field)
                if all_timesteps is not None:
                    self.nt = all_timesteps.shape[0]
                else:
                    # Fallback: count manually
                    if self.primary_field in self.store:
                        field_data = self.store[self.primary_field]
                        if hasattr(field_data, 'shape'):
                            self.nt = field_data.shape[0]
                        else:
                            self.nt = len(field_data.keys()) if hasattr(field_data, 'keys') else 1
                    else:
                        self.nt = 1
                
                # Get coordinates
                coords = self.data_manager.get_coordinates()
                self.x_coords = coords.get('x', np.linspace(-50, 50, self.nx))
                self.y_coords = coords.get('y', np.linspace(-50, 50, self.ny))
                self.z_coords = coords.get('z', np.linspace(0, 15, self.nz))
                
                print(f"Loaded data via DataManager: {self.nt} timesteps, {self.nx}x{self.ny}x{self.nz} grid")
                print(f"\n[FIELD INFO] Primary field: '{self.primary_field}'")
                print(f"  This is the main field being visualized. Secondary fields (qr, qv, u, v, w, etc.)")
                print(f"  will be combined with the primary field for richer visualization.")
                print(f"  To visualize a different field, use: --field <field_name>")
                print(f"  Available fields: {', '.join(self.data_manager.list_available_fields())}")
            except Exception as e:
                print(f"Warning: Could not use DataManager, falling back to direct zarr access: {e}")
                use_fallback = True
        else:
            use_fallback = True
        
        if use_fallback:
            # Fallback to original implementation
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
                    if len(z_coords_full) > self.nz:
                        self.z_coords = z_coords_full[:self.nz]
                    else:
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
        # Use DataManager if available for source-of-truth alignment
        if ENGINE_AVAILABLE and hasattr(self, 'data_manager'):
            try:
                # Load primary field
                all_primary = self.data_manager.get_all_timesteps(self.primary_field)
                if all_primary is not None:
                    self.primary_data = all_primary
                else:
                    # Fallback to direct access
                    primary_field_data = self.store[self.primary_field]
                    if hasattr(primary_field_data, 'shape'):
                        self.primary_data = primary_field_data[:]
                    else:
                        timestep_keys = sorted(primary_field_data.keys(), key=lambda x: int(x))
                        timestep_data = [primary_field_data[key][:] for key in timestep_keys]
                        self.primary_data = np.stack(timestep_data, axis=0)
                
                # Load secondary fields using DataManager's field discovery
                available_fields = self.data_manager.list_available_fields()
                self.secondary_fields = {}
                for field_name in available_fields:
                    if field_name != self.primary_field:
                        try:
                            field_all = self.data_manager.get_all_timesteps(field_name)
                            if field_all is not None:
                                self.secondary_fields[field_name] = field_all
                                print(f"Loaded field: {field_name} {field_all.shape}")
                        except Exception as e:
                            print(f"Warning: Could not load field {field_name}: {e}")
                            continue
                
                if self.secondary_fields:
                    print(f"Loaded secondary fields: {', '.join(self.secondary_fields.keys())}")
                else:
                    print("No secondary fields found - visualization will use primary field only")
                return
            except Exception as e:
                print(f"Warning: DataManager cache failed, using fallback: {e}")
        
        # Fallback to original implementation
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

        # Use ShaderManager if available
        if ENGINE_AVAILABLE:
            try:
                shader_dir = Path(__file__).parent / "shaders"
                self.shader_manager = ShaderManager(shader_dir, self.ctx)
                self.volume_prog = self.shader_manager.load_shader("volume", "volume_color")
                print("Loaded shaders via ShaderManager")
            except Exception as e:
                print(f"Warning: ShaderManager failed, using fallback: {e}")
                # Fallback to original
                shader_dir = Path(__file__).parent / "shaders"
                vertex_shader = self.load_shader_file(shader_dir / "volume" / "volume.vert")
                fragment_shader = self.load_shader_file(shader_dir / "volume" / "volume_color.frag")
                self.volume_prog = self.ctx.program(
                    vertex_shader=vertex_shader,
                    fragment_shader=fragment_shader
                )
        else:
            # Fallback to original implementation
            shader_dir = Path(__file__).parent / "shaders"
            # Try new location first, then fallback to old
            vert_path = shader_dir / "volume" / "volume.vert"
            frag_path = shader_dir / "volume" / "volume_color.frag"
            if not vert_path.exists():
                vert_path = shader_dir / "volume.vert"
            if not frag_path.exists():
                frag_path = shader_dir / "volume.frag"
            
            vertex_shader = self.load_shader_file(vert_path)
            fragment_shader = self.load_shader_file(frag_path)

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
            source = shader_path.read_text()
            # Resolve includes if using ShaderManager
            if ENGINE_AVAILABLE and hasattr(self, 'shader_manager'):
                try:
                    return self.shader_manager._resolve_includes(source, shader_path)
                except:
                    pass
            return source
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
        # DataManager returns (time, z, y, x), so [timestep] gives (z, y, x)
        primary = self.primary_data[timestep].copy()
        
        # Verify shape matches expectations
        if primary.shape != (self.nz, self.ny, self.nx):
            print(f"Warning: Primary data shape {primary.shape} doesn't match expected (nz={self.nz}, ny={self.ny}, nx={self.nx})")
            # Try to fix common issues
            if primary.shape == (self.nx, self.ny, self.nz):
                print("  Data appears to be (x, y, z), transposing to (z, y, x)...")
                primary = np.transpose(primary, (2, 1, 0))

        # Handle NaN and inf values
        nan_mask = np.isnan(primary)
        inf_mask = np.isinf(primary)
        invalid_mask = nan_mask | inf_mask
        
        if np.any(invalid_mask):
            invalid_count = np.sum(invalid_mask)
            print(f"Warning: timestep {timestep} has {invalid_count} invalid values (NaN/inf)")
            # Replace invalid values with a reasonable default based on field type
            if self.primary_field == 'theta':
                default_value = 300.0  # K
            elif self.primary_field in ['qr', 'qc', 'qv', 'qh', 'qg']:
                default_value = 0.0
            else:
                default_value = 0.0
            primary = np.where(invalid_mask, default_value, primary)

        # Get valid range (excluding any remaining invalid values)
        valid_data = primary[~np.isnan(primary) & ~np.isinf(primary)]
        if len(valid_data) == 0:
            print(f"Error: timestep {timestep} has no valid data")
            # Return empty RGBA
            return np.zeros((self.nx, self.ny, self.nz, 4), dtype=np.float32)
        
        primary_min, primary_max = valid_data.min(), valid_data.max()

        # Check for data corruption and apply field-specific normalization
        if self.primary_field == 'theta':
            # Theta should be reasonable temperature range (potential temperature)
            # Handle extreme negative values (might be data format issue)
            if primary_min < -1000 or primary_max > 1000 or not np.isfinite(primary_min) or not np.isfinite(primary_max):
                print(f"Warning: timestep {timestep} has suspicious theta range: {primary_min:.1f} to {primary_max:.1f}")
                
                # If values are extremely negative, might need to add offset or use absolute value
                if primary_min < -1000:
                    print(f"  Detected extreme negative values. Checking if data needs offset correction...")
                    # Try to find a reasonable offset (assume data should be around 300K)
                    # If all values are negative, add offset to bring them into reasonable range
                    if primary_max < 0:
                        offset = 300.0 - primary_max  # Shift so max becomes ~300
                        print(f"  Applying offset correction: +{offset:.1f} K")
                        primary = primary + offset
                        valid_data = primary[~np.isnan(primary) & ~np.isinf(primary)]
                        if len(valid_data) > 0:
                            primary_min, primary_max = valid_data.min(), valid_data.max()
                            print(f"  After offset: range [{primary_min:.1f}, {primary_max:.1f}] K")
                
                # Don't clamp yet - first check if we have a valid range for normalization
                # Use percentiles to handle outliers instead of hard clamping
                valid_data = primary[~np.isnan(primary) & ~np.isinf(primary)]
                if len(valid_data) > 0:
                    # Use percentiles to get robust min/max (ignore extreme outliers)
                    p1, p99 = np.percentile(valid_data, [1, 99])
                    if p99 - p1 > 1.0:  # If we have reasonable spread
                        # Use percentile range for normalization (more robust)
                        primary_min, primary_max = p1, p99
                        print(f"  Using percentile range (1-99%): [{primary_min:.1f}, {primary_max:.1f}] K")
                    else:
                        # If data is too compressed, try median ± std
                        median = np.median(valid_data)
                        std = np.std(valid_data)
                        if std > 0:
                            primary_min, primary_max = median - 3*std, median + 3*std
                            print(f"  Using median±3σ range: [{primary_min:.1f}, {primary_max:.1f}] K")
                        else:
                            # Last resort: use actual min/max
                            primary_min, primary_max = valid_data.min(), valid_data.max()
                            print(f"  Using actual range: [{primary_min:.1f}, {primary_max:.1f}] K")
                else:
                    primary_min, primary_max = 300.0, 350.0  # Fallback range
                    print(f"  Using fallback range: [{primary_min:.1f}, {primary_max:.1f}] K")

            # Standard linear normalization for temperature
            if primary_max > primary_min and np.isfinite(primary_min) and np.isfinite(primary_max):
                primary_norm = (primary - primary_min) / (primary_max - primary_min)
                primary_norm = np.clip(primary_norm, 0.0, 1.0)
            else:
                # Fallback: use reasonable default range
                print(f"Warning: Invalid range for theta, using default normalization")
                primary_norm = np.full_like(primary, 0.5, dtype=np.float32)

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
            # Default linear normalization for other fields (qr, qv, qc, etc.)
            # Handle case where field is all zeros (e.g., no rain at early timesteps)
            if primary_max == primary_min == 0.0:
                print(f"Warning: Field '{self.primary_field}' is all zeros at timestep {timestep}")
                print(f"  This is normal for fields like 'qr' (rain) at early timesteps.")
                print(f"  Visualization will use secondary fields (wind, theta, etc.) instead.")
                # Set to zero - secondary fields will provide the visualization
                primary_norm = np.zeros_like(primary)
            elif primary_max > primary_min and np.isfinite(primary_min) and np.isfinite(primary_max):
                primary_norm = (primary - primary_min) / (primary_max - primary_min)
                # Clamp to [0, 1] to handle any remaining edge cases
                primary_norm = np.clip(primary_norm, 0.0, 1.0)
            else:
                # Fallback: use a simple normalization if range is invalid
                if np.isfinite(primary_max) and np.isfinite(primary_min):
                    primary_norm = np.zeros_like(primary)
                else:
                    # If we still have inf/nan, use a default
                    primary_norm = np.full_like(primary, 0.5)
                    print(f"Warning: Using default normalization for timestep {timestep}")

        # Initialize RGBA channels
        r = primary_norm.copy()  # Red channel for primary field
        
        # If primary field is all zeros, that's okay - we'll use secondary fields
        # But add a small signal to red channel based on secondary fields so something shows
        if r.mean() == 0.0 and r.max() == 0.0:
            # Primary field is truly zero (e.g., no rain yet)
            # Use secondary fields to provide some red channel signal
            if 'theta' in self.secondary_fields:
                # Use theta as fallback for red channel
                theta_data = self.secondary_fields['theta'][timestep]
                theta_valid = theta_data[~np.isnan(theta_data) & ~np.isinf(theta_data)]
                if len(theta_valid) > 0:
                    theta_min, theta_max = theta_valid.min(), theta_valid.max()
                    if theta_max > theta_min:
                        theta_norm = np.clip((theta_data - theta_min) / (theta_max - theta_min), 0.0, 1.0)
                        r = np.maximum(r, theta_norm * 0.3)  # Use 30% of theta for red channel
                        print(f"  Using theta as fallback for red channel (primary field is zero)")
        
        g = np.zeros_like(primary_norm)  # Green for secondary fields
        b = np.zeros_like(primary_norm)  # Blue for wind/turbulence
        a = np.zeros_like(primary_norm)  # Alpha for opacity
        
        # Add primary field to alpha channel for visibility
        a = np.maximum(a, r * 0.3)  # Primary field contributes to opacity

        # Ensure we have valid RGBA data
        r = np.nan_to_num(r, nan=0.0, posinf=1.0, neginf=0.0)
        g = np.nan_to_num(g, nan=0.0, posinf=1.0, neginf=0.0)
        b = np.nan_to_num(b, nan=0.0, posinf=1.0, neginf=0.0)
        a = np.nan_to_num(a, nan=0.0, posinf=1.0, neginf=0.0)
        
        # Clamp all channels to [0, 1]
        r = np.clip(r, 0.0, 1.0)
        g = np.clip(g, 0.0, 1.0)
        b = np.clip(b, 0.0, 1.0)
        a = np.clip(a, 0.0, 1.0)
        
        # Debug: Check RGBA stats (only for first timestep)
        if timestep == self.valid_timesteps[0]:
            print(f"\n[TRANSFER FUNCTION] Timestep {timestep} ({self.primary_field}):")
            print(f"  Primary field range: [{primary_min:.2f}, {primary_max:.2f}]")
            
            # Warn if primary field is zero and suggest alternatives
            if primary_max == 0.0 and self.primary_field in ['qr', 'qc', 'qh', 'qg']:
                print(f"  ⚠️  Primary field '{self.primary_field}' is zero at this timestep (no precipitation yet)")
                print(f"  💡 Try: --field theta (thermal structure) or --field reflectivity_dbz (radar)")
                print(f"  💡 Or wait for later timesteps when precipitation forms")
            
            print(f"  RGBA stats: R:{r.mean():.3f}, G:{g.mean():.3f}, B:{b.mean():.3f}, A:{a.mean():.3f}")
            print(f"  RGBA ranges: R:[{r.min():.3f}, {r.max():.3f}], A:[{a.min():.3f}, {a.max():.3f}]")
            print(f"  Secondary fields used: {list(self.secondary_fields.keys())}")
            
            # Check which secondary fields actually contribute
            if 'qr' in self.secondary_fields:
                qr_data = self.secondary_fields['qr'][timestep]
                if qr_data.max() > 0:
                    print(f"    qr: range [{qr_data.min():.6f}, {qr_data.max():.6f}], mean={qr_data.mean():.6f}")
            if 'qv' in self.secondary_fields:
                qv_data = self.secondary_fields['qv'][timestep]
                print(f"    qv: range [{qv_data.min():.6f}, {qv_data.max():.6f}], mean={qv_data.mean():.6f}")
            
            # Calculate wind speed (handling case where primary might be a wind component)
            u_wind = primary if self.primary_field == 'u' else (self.secondary_fields.get('u', [None])[timestep] if 'u' in self.secondary_fields else None)
            v_wind = primary if self.primary_field == 'v' else (self.secondary_fields.get('v', [None])[timestep] if 'v' in self.secondary_fields else None)
            w_wind = primary if self.primary_field == 'w' else (self.secondary_fields.get('w', [None])[timestep] if 'w' in self.secondary_fields else None)
            
            if u_wind is not None and v_wind is not None and w_wind is not None:
                wind_speed = np.sqrt(u_wind**2 + v_wind**2 + w_wind**2)
                print(f"    wind: speed range [{wind_speed.min():.2f}, {wind_speed.max():.2f}] m/s")
                if wind_speed.max() > 1000:
                    print(f"    ⚠️  Very high wind speeds detected - data may need unit conversion")
            print()

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

        # Wind field visualization (blue channel)
        # Handle case where primary field might be a wind component
        u_data = None
        v_data = None
        w_data = None
        
        # Get wind components (may be primary or secondary)
        if self.primary_field == 'u':
            u_data = primary  # Use primary field
            v_data = self.secondary_fields.get('v', [None])[timestep] if 'v' in self.secondary_fields else None
            w_data = self.secondary_fields.get('w', [None])[timestep] if 'w' in self.secondary_fields else None
        elif self.primary_field == 'v':
            u_data = self.secondary_fields.get('u', [None])[timestep] if 'u' in self.secondary_fields else None
            v_data = primary  # Use primary field
            w_data = self.secondary_fields.get('w', [None])[timestep] if 'w' in self.secondary_fields else None
        elif self.primary_field == 'w':
            u_data = self.secondary_fields.get('u', [None])[timestep] if 'u' in self.secondary_fields else None
            v_data = self.secondary_fields.get('v', [None])[timestep] if 'v' in self.secondary_fields else None
            w_data = primary  # Use primary field
        else:
            # Primary is not a wind component, get all from secondary
            if 'u' in self.secondary_fields and 'v' in self.secondary_fields and 'w' in self.secondary_fields:
                u_data = self.secondary_fields['u'][timestep]
                v_data = self.secondary_fields['v'][timestep]
                w_data = self.secondary_fields['w'][timestep]
        
        # Calculate wind speed if we have all components
        if u_data is not None and v_data is not None and w_data is not None:
            wind_speed = np.sqrt(u_data**2 + v_data**2 + w_data**2)
            # Handle extreme wind speeds (data might be in wrong units)
            wind_max = np.nanpercentile(wind_speed, 99)  # Use 99th percentile to ignore outliers
            if wind_max > 1000:  # Unrealistic wind speeds
                print(f"  Warning: Very high wind speeds detected (max={wind_max:.1f} m/s). Data may need scaling.")
                wind_norm = np.clip(wind_speed / wind_max, 0.0, 1.0)  # Normalize to max
            else:
                wind_norm = np.clip(wind_speed / 50.0, 0.0, 1.0)  # Scale for typical storm winds (50 m/s)

            b = np.maximum(b, wind_norm)  # Blue for wind intensity

            # Add some wind contribution to opacity
            condensate_opacity = np.maximum(condensate_opacity, wind_norm * 0.2)

        # Add condensate opacity to final alpha
        a = np.maximum(a, condensate_opacity)

        # Ensure minimum opacity for structure visibility
        # Use a gradient based on primary field intensity
        min_opacity = 0.05 + r * 0.15  # Minimum 0.05, up to 0.2 based on primary field
        a = np.maximum(a, min_opacity)
        
        # Additional fallback: if alpha is still too low everywhere, boost it
        if a.max() < 0.1:
            print(f"Warning: Alpha channel very low (max={a.max():.3f}), boosting visibility")
            # Boost alpha based on any non-zero signal
            signal_mask = (r > 0.01) | (g > 0.01) | (b > 0.01)
            a = np.where(signal_mask, np.maximum(a, 0.2), a)
            # Also ensure primary field contributes to alpha
            a = np.maximum(a, r * 0.3)

        # Stack into RGBA
        rgba = np.stack([r, g, b, a], axis=-1)
        
        # Data is currently (z, y, x, 4) but texture expects (x, y, z, 4)
        # Transpose from (z, y, x, 4) to (x, y, z, 4)
        rgba = np.transpose(rgba, (2, 1, 0, 3))  # (z, y, x, 4) -> (x, y, z, 4)

        return rgba.astype(np.float32)

    def update_volume_texture(self, timestep: int):
        """Update 3D texture with current timestep data"""
        # Use precomputed data for smooth animation
        if timestep in self.precomputed_rgba:
            rgba_data = self.precomputed_rgba[timestep]
        else:
            # Fallback to on-demand computation if not precomputed
            rgba_data = self.create_transfer_function(timestep)
        
        # Verify data shape matches texture expectations
        expected_shape = (self.nx, self.ny, self.nz, 4)
        if rgba_data.shape != expected_shape:
            print(f"ERROR: Data shape mismatch! Expected {expected_shape}, got {rgba_data.shape}")
            print(f"  This will cause rendering issues. Data may need to be transposed.")
            # Try to fix common shape issues
            if rgba_data.shape[:3] == (self.nz, self.ny, self.nx):
                print(f"  Attempting transpose from (nz, ny, nx) to (nx, ny, nz)...")
                rgba_data = np.transpose(rgba_data, (2, 1, 0, 3))
        
        # Debug output: show what's being sent to GPU (only once, during first call)
        if not hasattr(self, '_texture_debug_logged'):
            r_mean = rgba_data[:, :, :, 0].mean()
            g_mean = rgba_data[:, :, :, 1].mean()
            b_mean = rgba_data[:, :, :, 2].mean()
            a_mean = rgba_data[:, :, :, 3].mean()
            a_max = rgba_data[:, :, :, 3].max()
            a_min = rgba_data[:, :, :, 3].min()
            non_zero_alpha = np.sum(rgba_data[:, :, :, 3] > 0.01)
            total_voxels = rgba_data.shape[0] * rgba_data.shape[1] * rgba_data.shape[2]
            
            print(f"\n[GPU TEXTURE] First timestep preview (timestep {timestep}):")
            print(f"  Shape: {rgba_data.shape} (expected: {expected_shape})")
            print(f"  RGBA means: R={r_mean:.3f}, G={g_mean:.3f}, B={b_mean:.3f}, A={a_mean:.3f}")
            print(f"  Alpha range: [{a_min:.3f}, {a_max:.3f}]")
            print(f"  Non-zero alpha voxels: {non_zero_alpha}/{total_voxels} ({100*non_zero_alpha/total_voxels:.1f}%)")
            
            if r_mean < 0.01:
                print(f"  WARNING: Red channel (primary field '{self.primary_field}') is nearly zero!")
                print(f"  This suggests normalization failed. Theta values may need different handling.")
            if a_max < 0.05:
                print(f"  WARNING: Alpha values very low! Rendering may be invisible.")
                print(f"  Try increasing opacity (UP arrow) or using preset 2 (press '2')")
            
            self._texture_debug_logged = True
        
        self.volume_texture.write(rgba_data.tobytes())

    def create_view_matrix(self) -> np.ndarray:
        """Create view matrix for orbiting camera"""
        # Use Camera system if available
        if ENGINE_AVAILABLE and hasattr(self, '_orbit_camera'):
            return self._orbit_camera.get_view_matrix()
        
        # Fallback to original implementation
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
        # Use Camera system if available
        if ENGINE_AVAILABLE and hasattr(self, '_orbit_camera'):
            return self._orbit_camera.get_projection_matrix(self.wnd.size[0], self.wnd.size[1])
        
        # Fallback to original implementation
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

        # Track previous frame to detect changes (only update texture when frame actually changes)
        if not hasattr(self, '_prev_render_frame'):
            self._prev_render_frame = -1
        
        # Update animation
        if self.is_playing and self.valid_timesteps:
            frame_index = int(time * 30 * self.animation_speed) % len(self.valid_timesteps)
            new_frame = self.valid_timesteps[frame_index]
            
            # Only update texture if frame changed (avoid spam)
            if new_frame != self.current_frame:
                self.current_frame = new_frame
                # Update texture with precomputed data (silently, no debug output)
                if hasattr(self, 'precomputed_rgba') and self.current_frame in self.precomputed_rgba:
                    rgba_data = self.precomputed_rgba[self.current_frame]
                    self.volume_texture.write(rgba_data.tobytes())
                else:
                    # Fallback: compute on demand (but suppress debug output)
                    rgba_data = self.create_transfer_function(self.current_frame)
                    self.volume_texture.write(rgba_data.tobytes())

            # Slowly orbit camera
            self.camera_angle += frame_time * 0.2
        elif self._prev_render_frame != self.current_frame:
            # Frame changed but not playing - still update texture (silently)
            if hasattr(self, 'precomputed_rgba') and self.current_frame in self.precomputed_rgba:
                rgba_data = self.precomputed_rgba[self.current_frame]
                self.volume_texture.write(rgba_data.tobytes())
        
        self._prev_render_frame = self.current_frame

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

        # Debug: Print render info occasionally (first frame only)
        if not hasattr(self, '_first_render_logged'):
            print(f"\n[RENDER] Starting rendering:")
            print(f"  Current frame: {self.current_frame}")
            print(f"  Opacity scale: {self.opacity_scale}")
            print(f"  Brightness: {self.brightness}")
            print(f"  Texture bound: {self.volume_texture}")
            print(f"  Shader program: {self.volume_prog}")
            self._first_render_logged = True

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
            # Number keys for presets (use ord() for character codes)
            elif key == ord('1'):
                self.opacity_scale = 0.1
                self.brightness = 1.0
                print("Preset 1: Default settings")
            elif key == ord('2'):
                self.opacity_scale = 0.5
                self.brightness = 2.0
                print("Preset 2: High visibility")
            elif key == ord('3'):
                self.opacity_scale = 1.0
                self.brightness = 3.0
                print("Preset 3: Maximum visibility")
            elif key == ord('4'):
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
        
        # Update camera system if available
        if ENGINE_AVAILABLE and hasattr(self, '_orbit_camera'):
            self._orbit_camera.set_angle(self.camera_angle)
            self._orbit_camera.set_height(self.camera_height)

def main():
    """Main entry point"""
    import sys
    # Check for --list-fields before opening window
    if '--list-fields' in sys.argv:
        # Handle field listing without opening window
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', required=True)
        parser.add_argument('--list-fields', action='store_true')
        args, _ = parser.parse_known_args()
        
        # Use DataManager to list fields without OpenGL
        try:
            import sys
            from pathlib import Path
            vis_path = Path(__file__).parent
            if str(vis_path) not in sys.path:
                sys.path.insert(0, str(vis_path))
            from core.data_manager import DataManager
            data_mgr = DataManager(args.input)
            fields = data_mgr.list_available_fields()
            print("\nAvailable fields in dataset:")
            print("=" * 50)
            for field in sorted(fields):
                print(f"  {field}")
            print("\n" + "=" * 50)
            print("Use --field FIELD_NAME to visualize a specific field as primary")
            print("Examples:")
            print("  --field theta      # Potential temperature (default)")
            print("  --field qr          # Rain water mixing ratio")
            print("  --field qv          # Water vapor mixing ratio")
            print("  --field reflectivity_dbz  # Radar reflectivity")
            print("\nNote: The renderer automatically combines ALL fields:")
            print("  - Primary field (--field) → Red channel")
            print("  - Condensate fields (qr, qc, etc.) → Green channel")
            print("  - Wind fields (u, v, w) → Blue channel")
            return 0
        except Exception as e:
            print(f"Error reading data: {e}")
            import traceback
            traceback.print_exc()
            return 1
    
    # Check dependencies first
    if not OPENGL_AVAILABLE:
        # Try to show help if requested, otherwise show error
        import sys
        if '--help' in sys.argv or '-h' in sys.argv:
            # Create a simple parser to show help
            import argparse
            parser = argparse.ArgumentParser(description="3D Volume Renderer for Supercell Simulations")
            SupercellVolumeRenderer.add_arguments(parser)
            parser.print_help()
            print("\n" + "="*60)
            print("NOTE: This script requires OpenGL dependencies:")
            print("  pip install moderngl moderngl-window zarr")
            print("\nFor 2D animation without OpenGL:")
            print("  python visualization/create_animation.py --input <data> --field theta --output animation.gif")
        else:
            print("ERROR: OpenGL dependencies not installed!")
            print("Required packages: moderngl, moderngl-window, zarr")
            print("Install with: pip install moderngl moderngl-window zarr")
            print("\nFor a simple 2D animation without OpenGL, use:")
            print("  python visualization/create_animation.py --input <data_path> --field theta --output animation.gif")
        return 1
    
    mglw.run_window_config(SupercellVolumeRenderer)

if __name__ == "__main__":
    import sys
    sys.exit(main() or 0)
