#version 330

in vec3 in_position;
out vec3 worldPos;

uniform mat4 mvpMatrix;

void main() {
    worldPos = in_position;
    gl_Position = mvpMatrix * vec4(in_position, 1.0);
}
