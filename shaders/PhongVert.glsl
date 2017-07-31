#version 410

/// model matrix
uniform mat4 m;
/// model view matrix
uniform mat4 mv;
/// model view projection matrix
uniform mat4 mvp;
/// normal matrix - from model to view space
uniform mat3 nm;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;

/// surface position in view space
out vec3 fr_pos;
/// surface normal in view space
out vec3 fr_normal;


void main()
{
	fr_pos = (mv * vec4(position, 1.0)).xyz;
	fr_normal = normalize(nm * normal);

	gl_Position = mvp * vec4(position, 1.0);
}
