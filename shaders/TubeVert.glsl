#version 410

layout ( location = 0 ) in vec3 position;
layout ( location = 1 ) in vec3 normal;

/// surface position in model space
out vec3 tcs_pos;
/// surface normal in model space
out vec3 tcs_norm;

void main()
{
	tcs_pos = position;
	tcs_norm = normalize( normal );

	gl_Position = vec4( 0.f );
}
