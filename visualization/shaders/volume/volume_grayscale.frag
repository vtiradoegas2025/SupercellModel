#version 330

in vec3 worldPos;
in vec3 viewDir;
out vec4 fragColor;

uniform sampler3D volumeTexture;
uniform vec3 cameraPos;
uniform float opacityScale;
uniform float brightness;
uniform float time;

// Include common utilities
#include "common/utils.glsl"

// Grayscale transfer function - maps scalar values to grayscale intensity
vec4 transferFunction(vec4 sample) {
    // Use red channel as primary intensity (or average if multiple channels)
    float intensity = sample.r;
    
    // Optional: use other channels for additional information
    // intensity = (sample.r + sample.g + sample.b) / 3.0;
    
    // Create grayscale color
    vec3 color = vec3(intensity, intensity, intensity);
    
    // Apply brightness
    color *= brightness;
    
    // Calculate opacity - higher intensity = more opaque
    float opacity = intensity * opacityScale;
    opacity = max(opacity, 0.05);  // Minimum visibility
    
    return vec4(color, opacity);
}

void main() {
    vec3 rayDir = normalize(viewDir);
    vec3 pos = cameraPos;
    
    vec4 accumulatedColor = vec4(0.0);
    float accumulatedAlpha = 0.0;
    
    // Maximum ray marching steps
    const int maxSteps = 200;
    const float stepSize = 0.015;
    
    // Ray march through volume
    for(int i = 0; i < maxSteps; i++) {
        // Transform world position to texture coordinates
        // Volume is assumed to be in [-1, 1] cube
        vec3 texCoord = (pos + 1.0) * 0.5;
        
        // Check if we're inside the volume
        if(all(greaterThanEqual(texCoord, vec3(0.0))) &&
           all(lessThanEqual(texCoord, vec3(1.0)))) {
            
            // Sample volume texture
            vec4 sampleColor = texture(volumeTexture, texCoord);
            
            // Apply transfer function
            vec4 finalColor = transferFunction(sampleColor);
            
            // Front-to-back blending
            float alpha = finalColor.a * (1.0 - accumulatedAlpha);
            accumulatedColor.rgb += finalColor.rgb * alpha;
            accumulatedAlpha += alpha;
            
            // Early termination
            if(accumulatedAlpha > 0.98) {
                break;
            }
        } else {
            // Outside volume - check if we should continue ray marching
            vec3 distToBox = max(abs(pos) - 1.0, 0.0);
            if(length(distToBox) > stepSize * 2.0) {
                break;
            }
        }
        
        // Step along ray
        pos += rayDir * stepSize;
    }
    
    // Ensure we have some minimum alpha for visibility
    accumulatedColor.a = accumulatedAlpha;
    
    // Gamma correction for display
    accumulatedColor.rgb = linearToSRGB(accumulatedColor.rgb);
    
    fragColor = accumulatedColor;
}
