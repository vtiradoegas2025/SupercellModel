#version 330

in vec3 worldPos;
in vec3 viewDir;
out vec4 fragColor;

uniform sampler3D volumeTexture;
uniform vec3 cameraPos;
uniform float opacityScale;
uniform float brightness;
uniform float time;  // For animation effects

// Transfer function - maps scalar values to colors
vec4 transferFunction(vec4 sample) {
    // RGBA channels encode different atmospheric fields:
    // sample.r = primary field (theta/buoyancy/temperature)
    // sample.g = condensate (qr + qc + cloud fields)
    // sample.b = water vapor (qv) or wind speed
    // sample.a = opacity (based on condensate + minimum visibility)

    vec3 color = vec3(0.0);

    // Primary field (red channel) - buoyancy/temperature
    // Red for warm updrafts, blue for cold downdrafts
    float temp = sample.r;
    color += mix(vec3(0.0, 0.0, 1.0), vec3(1.0, 0.0, 0.0), temp) * 0.6;

    // Condensate (green channel) - clouds and precipitation
    // White/grey for clouds and rain
    float condensate = sample.g;
    color += vec3(0.8, 0.8, 0.8) * condensate * 0.8;

    // Secondary field (blue channel) - water vapor or wind
    // Cyan for water vapor, magenta for wind
    float secondary = sample.b;
    color += vec3(0.0, 0.8, 1.0) * secondary * 0.4;

    // Apply brightness
    color *= brightness;

    // Calculate opacity - higher for condensate, minimum for visibility
    float opacity = sample.a * opacityScale;
    opacity = max(opacity, 0.02);

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
            // For a cube volume, we can break if we're completely outside
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

    // Gamma correction (approximate)
    accumulatedColor.rgb = pow(accumulatedColor.rgb, vec3(1.0/2.2));

    fragColor = accumulatedColor;
}
