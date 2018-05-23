#version 410

#define TWO_PI 6.283185307179586
#define ERR 1e-6f
#define DT 0.001

layout( quads, equal_spacing, cw ) in;

/// model view matrix
uniform mat4 mv;
/// model view projection matrix
uniform mat4 mvp;

/// tube radius
uniform float radius;

/// surface position in model space
in vec3 tes_pos[];
/// surface normal in model space
in vec3 tes_norm[];

/// surface position in view space
out vec3 fr_pos;
/// surface normal in view space
out vec3 fr_norm;
out vec2 fr_uv;

vec3 catmull_rom_cubic( float t )
{
	return 0.5f * ( ( 2.f * tes_pos[1 ])
			+ ( tes_pos[2] - tes_pos[0 ]) * t
			+ ( 2.f * tes_pos[0] - 5.f * tes_pos[1] + 4.f * tes_pos[2] - tes_pos[3] ) * t * t
			+ ( 3.f * tes_pos[1] - tes_pos[0] - 3.f * tes_pos[2] + tes_pos[3] ) * t * t * t );
}

vec3 catmull_rom_cubic_normal( float t )
{
	return 0.5f * ( ( 2.f * tes_norm[1] )
			+ ( tes_norm[2] - tes_norm[0] ) * t
			+ ( 2.f * tes_norm[0] - 5.f * tes_norm[1] + 4.f * tes_norm[2] - tes_norm[3] ) * t * t
			+ ( 3.f * tes_norm[1] - tes_norm[0] - 3.f * tes_norm[2] + tes_norm[3] ) * t * t * t );
}

// NOTE: n0, n1 need to be normalized vectors and 0 <= t <= 1
vec3 geometric_slerp( vec3 n0, vec3 n1, float t )
{
	float w = dot( n0, n1 );
	if ( abs( 1.f - w ) < ERR )
	{
		return n0;
	}

	w = acos( w );
	return normalize( ( sin( ( 1.f - t ) * w ) / sin( w ) ) * n0 + ( sin( t * w ) / sin( w ) ) * n1 );
}

vec3 slerp_normal( float t )
{
	float ei_1 = length( tes_pos[1] - tes_pos[0] );
	float ei = length( tes_pos[2] - tes_pos[1] );
	float ti = ei_1 / ( ei_1 + ei );
	vec3 nvi0 = geometric_slerp( tes_norm[0], tes_norm[1], ti );

	ei_1 = ei;
	ei = length( tes_pos[3] - tes_pos[2] );
	ti = ( ei_1 / ( ei_1 + ei ) );
	vec3 nvi1 = ( ( 1.f - ti ) < ERR )? tes_norm[1] : geometric_slerp( tes_norm[1], tes_norm[2], ti );

	return geometric_slerp( nvi0, nvi1, t );
}

void main ()
{
	float u = gl_TessCoord.x;
	float v = gl_TessCoord.y;
	float theta = u * TWO_PI;

	vec3 center = catmull_rom_cubic( v );
	vec3 tangent = ( catmull_rom_cubic( max( 0.f, v - 0.5f * DT ) )
					 - catmull_rom_cubic( min( v + 0.5f * DT, 1.f ) ) ) / DT;

//	vec3 normal = tes_norm[1];
//	vec3 normal = catmull_rom_cubic_normal( v );
	vec3 normal = slerp_normal( v );
	vec3 binormal = normalize( cross( tangent, normal ) );
	normal = cos( theta ) * normal + sin( theta ) * binormal;
	vec3 pos = center + normal * radius;

	fr_pos = ( mv * vec4( pos, 1.f ) ).xyz;
	fr_norm = ( transpose( inverse( mv ) ) * vec4( normal, 0.f ) ).xyz;
	fr_uv = gl_TessCoord.xy;

	gl_Position = mvp * vec4( pos, 1.f );
}
