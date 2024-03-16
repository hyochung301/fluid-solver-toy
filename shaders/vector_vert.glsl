#version 410 core

layout (location = 0) in vec2 aPos;

out vec2 iPos;

void main() {
    iPos = aPos;
    gl_Position = vec4(aPos, 1.0f, 1.0f);
}
