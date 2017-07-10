#version 410

layout (location = 0) in vec3 inVert;
layout (location = 1) in vec3 inKB;
layout (location = 2) in vec3 inM1;
layout (location = 3) in vec3 inM2;
layout (location = 4) in vec3 inForce;

out vec3 vert_gs;
out vec3 kb_gs;
out vec3 m1_gs;
out vec3 m2_gs;
out vec3 force_gs;

void main()
{
	vert_gs = inVert;
	kb_gs = inKB;
	m1_gs = inM1;
	m2_gs = inM2;
	force_gs = inForce;

	gl_Position = vec4(inVert, 1.0);
}
