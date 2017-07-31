#version 410

struct Light
{
	/// light position in view space
	vec3 position;
	/// light ambient emission color
	vec4 ambient;
	/// light diffuse emission color
	vec4 diffuse;
	/// light specular emission color
	vec4 specular;
};

struct Material
{
	/// ambient reflection color
	vec4 ambient;
	/// diffuse reflection color
	vec4 diffuse;
	/// specular reflection color
	vec4 specular;
	/// shininess - larger for surfaces that are smoother and more mirror-like
	float shininess;
};

/// point light
uniform Light light;
/// material
uniform Material mtl;

/// surface position in view space
in vec3 fr_pos;
/// surface normal in view space, normalized
in vec3 fr_normal;

/// output fragment colour
layout (location = 0) out vec4 fragColor;


void main ()
{
	// calculations take place in view space
	vec3 ldir = light.position - fr_pos;
	vec3 rdir = reflect(ldir, fr_normal);
	vec3 vdir = -fr_pos;

	fragColor = mtl.ambient * light.ambient +
			(mtl.diffuse * dot(ldir, fr_normal) * light.diffuse +
			 mtl.specular * pow(dot(rdir, vdir), mtl.shininess) * light.specular);
}
