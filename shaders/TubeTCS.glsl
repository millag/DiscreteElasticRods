#version 410

/// define the number of CPs in the output patch
layout ( vertices = 4 ) out;

/// attributes of the input CPs
in vec3 tcs_pos[];
in vec3 tcs_norm[];

/// attributes of the output CPs
out vec3 tes_pos[];
out vec3 tes_norm[];

void main ()
{
//	set tesselation levels only when first invocation of the shader
	if (gl_InvocationID == 0)
	{
		gl_TessLevelInner[0] = 6.f;
		gl_TessLevelOuter[1] = 6.f;
		gl_TessLevelOuter[3] = 6.f;

		gl_TessLevelInner[1] = 6.f;
		gl_TessLevelOuter[0] = 6.f;
		gl_TessLevelOuter[2] = 6.f;
	}

//	just copy input CPs to output
	tes_pos[gl_InvocationID] = tcs_pos[gl_InvocationID];
	tes_norm[gl_InvocationID] = tcs_norm[gl_InvocationID];
}
