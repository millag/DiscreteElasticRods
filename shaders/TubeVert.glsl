// glsl vertex shader
#version 410

layout (location = 0) in vec3 inVert;
layout (location = 1) in vec2 inUV;
layout (location = 2) in vec3 inNormal;

uniform mat4 M;

out vec3 wvert_cs;
out vec2 uv_cs;
out vec3 wnormal_cs;

void main()
{
    gl_Position = M * vec4(inVert, 1.0);
    wvert_cs = gl_Position.xyz;
    uv_cs = inUV;
    wnormal_cs = normalize(M * vec4(inNormal, 0.0)).xyz;
}
