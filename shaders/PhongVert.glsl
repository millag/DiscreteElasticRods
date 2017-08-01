#version 410

/// model view matrix
uniform mat4 mv;
/// model view projection matrix
uniform mat4 mvp;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;

/// surface position in view space
out vec3 fr_pos;
/// surface normal in view space
out vec3 fr_normal;


void main()
{
	fr_pos = (mv * vec4(position, 1.0)).xyz;
	fr_normal = (transpose(inverse(mv)) * vec4(normal, 0.0)).xyz;

	gl_Position = mvp * vec4(position, 1.0);
}
