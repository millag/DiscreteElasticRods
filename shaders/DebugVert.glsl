// glsl vertex shader
#version 410

layout (location = 0) in vec3 inVert;
layout (location = 1) in vec3 inKB;
layout (location = 2) in vec3 inM1;
layout (location = 3) in vec3 inM2;

out vec3 vert_gs;
out vec3 kb_gs;
out vec3 m1_gs;
out vec3 m2_gs;

void main()
{
    gl_Position = vec4(inVert, 1.0);
    vert_gs = gl_Position.xyz;
    kb_gs = inKB;
    m1_gs = inM1;
    m2_gs = inM2;
}
