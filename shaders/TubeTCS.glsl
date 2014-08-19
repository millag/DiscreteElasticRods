// glsl tessellation control shader
#version 410


// define the number of CPs in the output patch
layout (vertices = 4) out;

// attributes of the input CPs
in vec3 vert_cs[];
in vec3 normal_cs[];

// attributes of the output CPs
out vec3 vert_es[];
out vec3 normal_es[];

void main ()
{
//     set tesselation levels only when first invocation of the shader
    if (gl_InvocationID == 0)
    {
        gl_TessLevelInner[0] = 6.0;
        gl_TessLevelOuter[1] = 6.0;
        gl_TessLevelOuter[3] = 6.0;

        gl_TessLevelInner[1] = 5.0;
        gl_TessLevelOuter[0] = 5.0;
        gl_TessLevelOuter[2] = 5.0;
    }

//     Just pass the control points of the output patch
    vert_es[gl_InvocationID] = vert_cs[gl_InvocationID];
    normal_es[gl_InvocationID] = normal_cs[gl_InvocationID];
}
