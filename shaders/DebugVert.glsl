#version 410

layout ( location = 0 ) in vec3 position;
layout ( location = 1 ) in vec3 kb;
layout ( location = 2 ) in vec3 m1;
layout ( location = 3 ) in vec3 m2;
layout ( location = 4 ) in vec3 force;

out vec3 gs_pos;
out vec3 gs_kb;
out vec3 gs_m1;
out vec3 gs_m2;
out vec3 gs_force;

void main()
{
	gs_pos = position;
	gs_kb = kb;
	gs_m1 = m1;
	gs_m2 = m2;
	gs_force = force;

	gl_Position = vec4( 0.f );
}
