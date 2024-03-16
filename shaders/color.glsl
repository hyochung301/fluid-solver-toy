#version 410 core

out vec4 outColor;

in vec2 iPos;

void main(){
    vec3 c = vec3( 1.f, 0.f, 0.f);
    outColor = vec4(c,0.);
}
