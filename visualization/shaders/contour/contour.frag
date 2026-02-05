#version 330

in vec3 worldPos;
out vec4 fragColor;

uniform sampler3D volumeTexture;
uniform float contourSpacing;
uniform float contourWidth;
uniform vec3 contourColor;
uniform float opacity;

// Dense contour rendering for turbulence visualization
// Similar to Lewellen et al. (2008) style

void main() 
{
    // Transform world position to texture coordinates
    vec3 texCoord = (worldPos + 1.0) * 0.5;
    
    // Check if we're inside the volume
    if(any(lessThan(texCoord, vec3(0.0))) || any(greaterThan(texCoord, vec3(1.0)))) 
    {
        fragColor = vec4(0.0);
        return;
    }
    
    // Sample volume texture
    float value = texture(volumeTexture, texCoord).r;
    
    // Calculate contour lines
    // Use modulo to create evenly spaced contours
    float contourValue = mod(value, contourSpacing);
    
    // Create contour lines with smooth edges
    float line = smoothstep(contourWidth * 0.5, contourWidth, abs(contourValue - contourSpacing * 0.5));
    
    // Invert so lines are visible (1.0 = line, 0.0 = background)
    line = 1.0 - line;
    
    // Apply contour color
    vec3 color = contourColor * line;
    
    fragColor = vec4(color, opacity * line);
}
