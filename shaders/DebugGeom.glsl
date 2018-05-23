#version 410

struct Light
{
/// light ambient emission color
	vec3 ambient;
/// light diffuse emission color
	vec3 diffuse;
/// light specular emission color
	vec3 specular;
/// light position in view space
	vec3 position;
};

/// point light
uniform Light light;
/// model view matrix
uniform mat4 mv;
/// model view projection matrix
uniform mat4 mvp;


layout ( lines_adjacency ) in;
layout ( line_strip, max_vertices = 10 ) out;

in vec3 gs_pos[];
in vec3 gs_kb[];
in vec3 gs_m1[];
in vec3 gs_m2[];
in vec3 gs_force[];

out vec3 edir;
out vec3 ldir;
out vec3 tdir;
out vec4 color;

void main ()
{
	mat4 imv = inverse( mv );
	vec3 eyePos = ( imv * vec4( 0.f, 0.f, 0.f, 1.f ) ).xyz;
	vec3 lightPos = ( imv * vec4( light.position, 1.f ) ).xyz;

//	draw line
	edir = normalize( eyePos - gs_pos[1] );
	ldir = normalize( lightPos - gs_pos[1] );
	tdir = normalize( gs_pos[2] - gs_pos[1] );
	color = vec4( 1.f, 1.f, 0.f, 1.f );
	gl_Position = mvp * vec4( gs_pos[1], 1.f );
	EmitVertex();

	edir = normalize( eyePos - gs_pos[2] );
	ldir = normalize( lightPos - gs_pos[2] );
	color = vec4( 1.f, 1.f, 0.f, 1.f );
	gl_Position = mvp * vec4( gs_pos[2], 1.f );
	EmitVertex();

	EndPrimitive();

//	draw kb
	edir = normalize( eyePos - gs_pos[1] );
	ldir = normalize( lightPos - gs_pos[1] );
	tdir = normalize( gs_kb[1] );
	color = vec4( 0.f, 1.f, 0.f, 1.f );
	gl_Position = mvp * vec4( gs_pos[1], 1.f );
	EmitVertex();

	vec3 pos1 = gs_pos[1] + tdir;
	edir = normalize( eyePos - pos1 );
	ldir = normalize( lightPos - pos1 );
	color = vec4( 0.f, 1.f, 0.f, 1.f );
	gl_Position = mvp * vec4( pos1, 1.f );
	EmitVertex();

	EndPrimitive();

//	draw force
	edir = normalize( eyePos - gs_pos[1] );
	ldir = normalize( lightPos - gs_pos[1] );
	tdir = normalize( gs_force[1] );
	color = vec4( 0.f, 1.f, 1.f, 1.f );
	gl_Position = mvp * vec4( gs_pos[1], 1.f );
	EmitVertex();

	pos1 = gs_pos[1] + tdir;
	edir = normalize( eyePos - pos1 );
	ldir = normalize( lightPos - pos1 );
	color = vec4( 0.f, 1.f, 1.f, 1.f );
	gl_Position = mvp * vec4( pos1, 1.f );
	EmitVertex();

	EndPrimitive();

//  draw material frame m1, m2
	pos1 = ( gs_pos[1] + gs_pos[2] ) * 0.5f;
	const float len = 0.1f;

	edir = normalize( eyePos - pos1 );
	ldir = normalize( lightPos - pos1 );
	tdir = normalize( gs_m1[1] );
	color = vec4( 1.f, 0.f, 0.f, 1.f);
	gl_Position = mvp * vec4( pos1, 1.f );
	EmitVertex();

	vec3 pos2 = pos1 + len * tdir;
	edir = normalize( eyePos - pos2 );
	ldir = normalize( lightPos - pos2 );
	color = vec4( 1.f, 0.f, 0.f, 1.f );
	gl_Position = mvp * vec4( pos2, 1.f );
	EmitVertex();

	EndPrimitive();

	edir = normalize( eyePos - pos1 );
	ldir = normalize( lightPos - pos1 );
	tdir = normalize( gs_m2[1] );
	color = vec4( 0.f, 0.f, 1.f, 1.f );
	gl_Position = mvp * vec4( pos1, 1.f );
	EmitVertex();

	pos2 = pos1 + len * tdir;
	edir = normalize( eyePos - pos2 );
	ldir = normalize( lightPos - pos2 );
	color = vec4( 0.f, 0.f, 1.f, 1.f );
	gl_Position = mvp * vec4( pos2, 1.f );
	EmitVertex();

	EndPrimitive();
}
