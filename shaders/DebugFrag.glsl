#version 410

layout (location = 0) out vec4 fragColour;

in vec3 eDir;
in vec3 lDir;
in vec3 tDir;
in vec4 color;


void main()
{
	vec3 e = normalize(eDir);
	vec3 l = normalize(lDir);
	vec3 t = normalize(tDir);

	float dotTL = dot(l, t);
	float dotTE = dot(t, e);
	vec4 diffuse = color * sqrt(1.0 - dotTL * dotTL);
	vec4 specular = vec4(1.,1.,1.,1.) *  pow((dotTL * dotTE + sqrt(1.0 - dotTL * dotTL) * sqrt(1.0 - dotTE * dotTE)), 0.2);
	fragColour = diffuse;
}

