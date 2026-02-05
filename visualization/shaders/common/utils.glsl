// Common utility functions for shaders

// Convert RGB to grayscale
float rgbToGrayscale(vec3 rgb) 
{
    return dot(rgb, vec3(0.299, 0.587, 0.114));
}

// Smooth step with configurable edge
float smoothStep(float edge0, float edge1, float x) 
{
    float t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

// Gamma correction
vec3 gammaCorrect(vec3 color, float gamma) 
{
    return pow(color, vec3(1.0 / gamma));
}

// Linear to sRGB
vec3 linearToSRGB(vec3 linear) 
{
    return gammaCorrect(linear, 2.2);
}
