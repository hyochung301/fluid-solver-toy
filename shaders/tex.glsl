#version 410 core

uniform sampler2D uTexslot;
out vec4 outColor;
in vec2 iUV;

void main(){
    vec4 c = texture(uTexslot, iUV);
    outColor = c;
}
