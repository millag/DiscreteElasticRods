// glsl tessellation evaluation shader
#version 410

#define TWO_PI 6.283185307179586
#define ERR 1e-6f
#define DT 0.001

layout(quads, equal_spacing, ccw) in;

uniform float Radius;
uniform mat4 M;
uniform mat4 MV;
uniform mat4 MVP;
uniform mat3 normalMatrix;

in vec3 vert_es[];
in vec3 normal_es[];

out vec3 vert_fr;
out vec3 tangent_fr;
out vec3 normal_fr;
out vec2 uv_fr;


vec3 catmull_rom_cubic(float t)
{
    return 0.5f * ((2 * vert_es[1]) +
                  (vert_es[2] - vert_es[0]) * t +
                  (2 * vert_es[0] - 5 * vert_es[1] + 4 * vert_es[2] - vert_es[3]) * t * t +
                  (3 * vert_es[1] - vert_es[0] - 3 * vert_es[2] + vert_es[3]) * t * t * t);
}

vec3 catmull_rom_cubic_normal(float t)
{
    return 0.5f * ((2 * normal_es[1]) +
                  (normal_es[2] - normal_es[0]) * t +
                  (2 * normal_es[0] - 5 * normal_es[1] + 4 * normal_es[2] - normal_es[3]) * t * t +
                  (3 * normal_es[1] - normal_es[0] - 3 * normal_es[2] + normal_es[3]) * t * t * t);
}

// NOTE: n0, n1 need to be normalized vectors and 0 <= t <= 1
vec3 geometric_slerp(vec3 n0, vec3 n1, float t)
{
    float w = dot(n0, n1);
    if (abs(1 - w) < ERR)
        return n0;

    w = acos(w);
    return normalize((sin((1 - t)*w) / sin(w)) * n0 + (sin(t * w) / sin(w)) * n1);
}

vec3 slerp_normal(float t)
{
    float ei_1 = length(vert_es[1] - vert_es[0]);
    float ei = length(vert_es[2] - vert_es[1]);
    float ti = ei_1 / (ei_1 + ei);
    vec3 nvi0 = geometric_slerp(normal_es[0], normal_es[1], ti);

    ei_1 = ei;
    ei = length(vert_es[3] - vert_es[2]);
    ti = (ei_1 / (ei_1 + ei));
    vec3 nvi1 = ((1 - ti) < ERR)? normal_es[1] : geometric_slerp(normal_es[1], normal_es[2], ti);

    return geometric_slerp(nvi0, nvi1, t);
}



void main ()
{
    float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;
    uv_fr = gl_TessCoord.xy;

    vec3 center = catmull_rom_cubic(v);
    tangent_fr = (catmull_rom_cubic( max(0.0, v - 0.5 * DT) ) - catmull_rom_cubic( min(v + 0.5 * DT, 1.0) )) / DT;

    vec3 normal = normal_es[1];
//    vec3 normal = catmull_rom_cubic_normal(v);
//    vec3 normal = slerp_normal(v);
    vec3 binormal = normalize( cross(tangent_fr, normal) );

    float theta = u * TWO_PI;
    normal_fr = cos(theta) * normal + sin(theta) * binormal;
    vert_fr = normal_fr * Radius + center;

    gl_Position = MVP * vec4(vert_fr, 1.0);
}
