#version 330

in vec3 in_position;
out vec3 worldPos;
out vec3 viewDir;

uniform mat4 mvpMatrix;
uniform vec3 cameraPos;

void main() {
    // Pass world position to fragment shader
    worldPos = in_position;

    // Calculate view direction from camera to vertex
    viewDir = normalize(in_position - cameraPos);

    // Transform to clip space
    gl_Position = mvpMatrix * vec4(in_position, 1.0);
}
