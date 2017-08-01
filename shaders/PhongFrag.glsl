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
/// surface normal in view space, normalized
in vec3 fr_normal;

/// output fragment colour
layout (location = 0) out vec4 fragColor;


void main ()
{
	// calculations take place in view space
	vec3 ldir = normalize(light.position - fr_pos);
	vec3 ndir = normalize(fr_normal);
	vec3 rdir = reflect(ldir, ndir);
	vec3 vdir = -normalize(fr_pos);

	vec3 color = light.ambient * mtl.ambient +
			(light.diffuse * mtl.diffuse * max(dot(ldir, ndir), 0.) +
			 light.specular * mtl.specular * pow(max(dot(rdir, vdir), 0.), mtl.shininess));

	fragColor = vec4(color, 1.0);
}
