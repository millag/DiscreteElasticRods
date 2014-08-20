//glsl franment shader
#version 410

layout (location = 0) out vec4 frag_colour;

struct Materials
{
  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  float shininess;
};

struct Lights
{
  vec4 position;
  vec4 ambient;
  vec4 diffuse;
  vec4 specular;
  float constantAttenuation;
  float spotCosCutoff;
  float quadraticAttenuation;
  float linearAttenuation;
};


uniform vec4 Colour;

uniform mat4 M;
uniform mat4 MV;
uniform mat4 MVP;
uniform mat3 normalMatrix;

uniform Lights light;
uniform Materials material;

in vec3 vert_fr;
in vec3 tangent_fr;
in vec3 normal_fr;
in vec2 uv_fr;


void main()
{
    vec3 L = light.position.xyz;
    vec3 N = normalize( normalMatrix * normal_fr );
    N = faceforward(N, vec3(0), N);
    vec3 V = (MV * vec4(vert_fr, 1.0)).xyz;
    vec3 H = normalize( L + V );

    vec4 ambient = material.ambient * light.ambient;
    vec4 diffuse = material.diffuse * light.diffuse * dot(N, L);
    vec4 specular = material.specular * light.specular * pow( max(dot(N, H), 0.0),  material.shininess);

    frag_colour = ambient + diffuse;// + specular;
}

