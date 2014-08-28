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


uniform vec4 Col1;
uniform vec4 Col2;

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
    vec3 N = -normalize( normalMatrix * normal_fr );
    vec3 V = (MV * vec4(vert_fr, 1.0)).xyz;
    vec3 H = normalize( L + V );


    float repeatCount = 15;
    float whichStripe = floor(uv_fr.x * repeatCount);
    vec4 col = mix(Col1, Col2, mod( whichStripe, 2));

    vec4 ambient = material.ambient * light.ambient;
    vec4 diffuse = material.diffuse * light.diffuse * dot(N, L) * col;
    vec4 specular = material.specular * light.specular * pow( max(dot(N, H), 0.0),  material.shininess);

    frag_colour = ambient + diffuse + specular;
}

