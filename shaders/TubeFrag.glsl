//glsl franment shader
#version 410

layout (location =0) out vec4 frag_colour;

uniform vec4 Colour;
uniform mat4 M;
uniform mat4 MV;
uniform mat4 MVP;
uniform mat3 normalMatrix;

in vec3 vert_fr;
in vec3 tangent_fr;
in vec3 normal_fr;
in vec2 uv_fr;


void main()
{
//    vec3 e = normalize(eDir);
//    vec3 l = normalize(lDir);
//    vec3 t = normalize(tDir);

//    float dotTL = dot(l, t);
//    float dotTE = dot(t, e);
//    vec4 diffuse = Colour * sqrt(1.0 - dotTL * dotTL);
//    vec4 specular = vec4(1.,1.,1.,1.) *  pow((dotTL * dotTE + sqrt(1.0 - dotTL * dotTL) * sqrt(1.0 - dotTE * dotTE)), 0.2);
//    frag_colour = diffuse;

    frag_colour = Colour;
}

