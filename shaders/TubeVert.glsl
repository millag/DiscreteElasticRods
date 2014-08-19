// glsl vertex shader
#version 410

layout (location = 0) in vec3 inVert;
layout (location = 1) in vec3 inNormal;

out vec3 vert_cs;
out vec3 normal_cs;

void main()
{
    vert_cs = inVert;
    normal_cs = normalize(inNormal);
}
