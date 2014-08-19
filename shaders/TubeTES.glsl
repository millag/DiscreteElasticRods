// glsl tessellation evaluation shader
#version 410

#define TWO_PI 6.283185307179586

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
//TODO: need look up for formula of quadratic catmull-rom
vec3 catmull_rom_quadratic(float t)
{
    return 0.5f * ((2 * normal_es[1]) +
                  (normal_es[2] - normal_es[0]) * t +
                  (2 * normal_es[0] - 5 * normal_es[1] + 4 * normal_es[2] - normal_es[3]) * t * t +
                  (3 * normal_es[1] - normal_es[0] - 3 * normal_es[2] + normal_es[3]) * t * t * t);
}

void main ()
{
    float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;
    uv_fr = gl_TessCoord.xy;

    vec3 center = catmull_rom_cubic(v);

    float dt = 0.001;
    tangent_fr = (catmull_rom_cubic( max(0.0, v - 0.5 * dt) ) - catmull_rom_cubic( min(v + 0.5 * dt, 1.0) )) / dt;
    normal_fr = catmull_rom_quadratic(v);

    vec3 tangent = normalize(tangent_fr);
    vec3 normal = normalize(normal_fr);
    vec3 binormal = normalize( cross(tangent, normal) );


    float theta = u * TWO_PI;
    vert_fr = (cos(theta) * normal + sin(theta) * binormal) * Radius + center;

    gl_Position = MVP * vec4(vert_fr, 1.0);
}
