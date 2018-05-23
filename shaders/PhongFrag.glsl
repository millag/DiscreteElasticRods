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

struct Material
{
/// ambient reflection color
	vec3 ambient;
/// diffuse reflection color
	vec3 diffuse;
/// specular reflection color
	vec3 specular;
/// shininess - larger for surfaces that are smoother and more mirror-like
	float shininess;
};

/// point light
uniform Light light;
/// material
uniform Material mtl;

/// surface position in view space
in vec3 fr_pos;
/// surface normal in view space
in vec3 fr_normal;

/// output fragment color
layout ( location = 0 ) out vec4 fragColor;

void main ()
{
// calculations are in view space
	vec3 ldir = normalize( light.position - fr_pos );
	vec3 ndir = normalize( fr_normal );
	vec3 vdir = normalize( -fr_pos );
	vec3 hdir = normalize( ldir + vdir );
//	vec3 rdir = reflect( -ldir, ndir );

	vec3 finalCol = light.ambient * mtl.ambient
			+ light.diffuse * mtl.diffuse * max( dot( ldir, ndir ), 0.f )
			+ light.specular * mtl.specular * pow( max( dot( hdir, ndir ), 0.f ), mtl.shininess );
//			+ light.specular * mtl.specular * pow( max( dot( rdir, vdir ), 0.f ), mtl.shininess / 4 );

	fragColor = vec4( finalCol, 1.f );
}
