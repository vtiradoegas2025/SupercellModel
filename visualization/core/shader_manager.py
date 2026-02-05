"""
Shader Manager for loading, compiling, and composing GLSL shaders.

Supports:
- Simple shader file loading
- #include directives for shared code
- Shader composition for advanced effects
"""

import moderngl as mgl
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import re


class ShaderManager:
    """
    Manages shader loading, compilation, and composition.
    
    Supports both simple shader swapping and advanced composition.
    """
    
    def __init__(self, shader_dir: Union[str, Path], ctx: Optional[mgl.Context] = None):
        """
        Initialize ShaderManager.
        
        Args:
            shader_dir: Directory containing shader files
            ctx: ModernGL context (if None, will be set when loading shaders)
        """
        self.shader_dir = Path(shader_dir)
        if not self.shader_dir.exists():
            raise FileNotFoundError(f"Shader directory does not exist: {self.shader_dir}")
        
        self.ctx = ctx
        self._shader_cache: Dict[str, str] = {}
        self._program_cache: Dict[str, mgl.Program] = {}
    
    def set_context(self, ctx: mgl.Context):
        """Set OpenGL context (required before loading shaders)"""
        self.ctx = ctx
    
    def _resolve_includes(self, source: str, base_path: Path) -> str:
        """
        Resolve #include directives in shader source.
        
        Supports:
        - #include "path/to/file.glsl"
        - #include <common/file.glsl>  (searches in common/ directory)
        
        Args:
            source: Shader source code
            base_path: Base path for resolving relative includes
            
        Returns:
            Source code with includes resolved
        """
        lines = source.split('\n')
        result_lines = []
        included_files = set()  # Prevent circular includes
        
        for line in lines:
            # Match #include directives
            include_match = re.match(r'#include\s+["<](.+?)[">]', line)
            if include_match:
                include_path = include_match.group(1)
                
                # Resolve path
                if include_path.startswith('common/'):
                    # Search in common directory
                    include_file = self.shader_dir / include_path
                else:
                    # Relative to base path
                    include_file = base_path.parent / include_path
                
                if include_file in included_files:
                    # Already included, skip
                    continue
                
                if not include_file.exists():
                    raise FileNotFoundError(f"Included shader file not found: {include_file}")
                
                # Read and recursively resolve includes
                included_content = include_file.read_text()
                included_files.add(include_file)
                resolved_content = self._resolve_includes(included_content, include_file)
                result_lines.append(f"// Included from {include_file.name}")
                result_lines.extend(resolved_content.split('\n'))
                result_lines.append(f"// End include {include_file.name}")
            else:
                result_lines.append(line)
        
        return '\n'.join(result_lines)
    
    def load_shader_source(self, shader_name: str, shader_type: str = 'frag') -> str:
        """
        Load shader source code with include resolution.
        
        Args:
            shader_name: Name of shader (without extension)
            shader_type: 'vert' or 'frag'
            
        Returns:
            Resolved shader source code
        """
        # Try different possible locations
        possible_paths = [
            self.shader_dir / f"{shader_name}.{shader_type}",
            self.shader_dir / shader_type / f"{shader_name}.{shader_type}",
            self.shader_dir / "volume" / f"{shader_name}.{shader_type}",
            self.shader_dir / "contour" / f"{shader_name}.{shader_type}",
        ]
        
        shader_path = None
        for path in possible_paths:
            if path.exists():
                shader_path = path
                break
        
        if shader_path is None:
            raise FileNotFoundError(
                f"Shader not found: {shader_name}.{shader_type}. "
                f"Searched in: {[str(p) for p in possible_paths]}"
            )
        
        # Check cache
        cache_key = str(shader_path)
        if cache_key in self._shader_cache:
            return self._shader_cache[cache_key]
        
        # Load and resolve includes
        source = shader_path.read_text()
        resolved_source = self._resolve_includes(source, shader_path)
        
        # Cache
        self._shader_cache[cache_key] = resolved_source
        
        return resolved_source
    
    def load_shader(self, vertex_shader: str, fragment_shader: str) -> mgl.Program:
        """
        Load and compile a shader program.
        
        Args:
            vertex_shader: Name of vertex shader (without extension)
            fragment_shader: Name of fragment shader (without extension)
            
        Returns:
            Compiled shader program
        """
        if self.ctx is None:
            raise RuntimeError("OpenGL context not set. Call set_context() first.")
        
        # Check cache
        cache_key = f"{vertex_shader}:{fragment_shader}"
        if cache_key in self._program_cache:
            return self._program_cache[cache_key]
        
        # Load shader sources
        vert_source = self.load_shader_source(vertex_shader, 'vert')
        frag_source = self.load_shader_source(fragment_shader, 'frag')
        
        # Compile program
        try:
            program = self.ctx.program(
                vertex_shader=vert_source,
                fragment_shader=frag_source
            )
        except Exception as e:
            raise RuntimeError(
                f"Failed to compile shader program ({vertex_shader}, {fragment_shader}): {e}"
            )
        
        # Cache
        self._program_cache[cache_key] = program
        
        return program
    
    def compose_shaders(self, shader_parts: List[str], shader_type: str = 'frag') -> str:
        """
        Compose multiple shader parts into a single shader.
        
        This is for advanced users who want to combine multiple shader effects.
        
        Args:
            shader_parts: List of shader names or file paths to compose
            shader_type: 'vert' or 'frag'
            
        Returns:
            Composed shader source code
        """
        composed_parts = []
        
        for part in shader_parts:
            if Path(part).exists():
                # Direct file path
                source = Path(part).read_text()
            else:
                # Shader name
                source = self.load_shader_source(part, shader_type)
            
            composed_parts.append(source)
        
        # Combine with separator
        separator = f"\n// === Composed from {len(shader_parts)} parts ===\n"
        return separator.join(composed_parts)
    
    def get_program(self, vertex_shader: str, fragment_shader: str) -> mgl.Program:
        """Get cached or load shader program"""
        return self.load_shader(vertex_shader, fragment_shader)
    
    def clear_cache(self):
        """Clear shader and program caches"""
        self._shader_cache.clear()
        self._program_cache.clear()
    
    def list_available_shaders(self) -> Dict[str, List[str]]:
        """
        List available shaders in the shader directory.
        
        Returns:
            Dictionary with 'vert' and 'frag' keys containing lists of shader names
        """
        shaders = {'vert': [], 'frag': []}
        
        # Search recursively
        for shader_file in self.shader_dir.rglob("*.vert"):
            rel_path = shader_file.relative_to(self.shader_dir)
            shaders['vert'].append(str(rel_path.with_suffix('')))
        
        for shader_file in self.shader_dir.rglob("*.frag"):
            rel_path = shader_file.relative_to(self.shader_dir)
            shaders['frag'].append(str(rel_path.with_suffix('')))
        
        return shaders
