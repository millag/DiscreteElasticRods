// glsl tessellation control shader
#version 410
#ension GL_ARB_tessellation_shader: enable


// define the number of CPs in the output patch
layout (vertices = 4) out;

// attributes of the input CPs
in vec3 wPos_cs[];
in vec2 uv_cs[];
in vec3 wNormal_cs[];

// attributes of the output CPs
out vec3 wPos_es[];
out vec2 uv_es[];
out vec3 wNormal_es[];

void main ()
{
    if (gl_InvocationID == 0)
    {
        gl_TessLevelOuter[0] = 1.0;
        gl_TessLevelOuter[1] = 5.0;
    }

    // Set the control points of the output patch
    wPos_es[gl_InvocationID] = wPos_cs[gl_InvocationID];
    uv_es[gl_InvocationID] = uv_cs[gl_InvocationID];
    wNormal_es[gl_InvocationID] = wNormal_cs[gl_InvocationID];
}
