#version 410

layout ( location = 0 ) out vec4 fragColour;

in vec3 edir;
in vec3 ldir;
in vec3 tdir;
in vec4 color;

void main()
{
	vec3 e = normalize( edir );
	vec3 l = normalize( ldir );
	vec3 t = normalize( tdir );

	float dotTL = dot( l, t );
	float dotTE = dot( t, e );
	vec4 diffuse = color * sqrt( 1.f - dotTL * dotTL );
	vec4 specular = vec4( 1.f, 1.f, 1.f, 1.f )
			* pow( ( dotTL * dotTE + sqrt( 1.f - dotTL * dotTL ) * sqrt( 1.f - dotTE * dotTE ) ), 0.2f );

	fragColour = diffuse;
}
