// GLSL GEOMETRY SHADER
#version 410

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

uniform mat4 MVP;
// array of lights
uniform Lights light;
uniform vec3 eyePos;
uniform vec4 Colour;

layout (lines_adjacency) in;
layout (line_strip, max_vertices = 10) out;

in vec3 vert_gs[];
in vec3 kb_gs[];
in vec3 m1_gs[];
in vec3 m2_gs[];
in vec3 force_gs[];

out vec3 eDir;
out vec3 lDir;
out vec3 tDir;
out vec4 color;

void main ()
{
//  draw line
    gl_Position = vec4(vert_gs[1], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    tDir = normalize(vert_gs[2] - vert_gs[1]);
    color = Colour;
    gl_Position = MVP * gl_Position;
    EmitVertex();

    gl_Position = vec4(vert_gs[2], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    tDir = normalize(vert_gs[2] - vert_gs[1]);
    color = Colour;
    gl_Position = MVP * gl_Position;
    EmitVertex();

    EndPrimitive();

    vec3 pos = (vert_gs[1] + vert_gs[2]) * 0.5;
    const float len = 0.1;
//  draw kb
    tDir = normalize(kb_gs[1]);
    gl_Position = vec4(vert_gs[1], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 1.0, 0.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    tDir = normalize(kb_gs[1]);
    gl_Position = vec4(vert_gs[1] + kb_gs[1], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 1.0, 0.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    EndPrimitive();

//  draw force
    tDir = normalize(force_gs[1]);
    gl_Position = vec4(vert_gs[1], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 1.0, 1.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    tDir = normalize(force_gs[1]);
    gl_Position = vec4(vert_gs[1] + force_gs[1], 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 1.0, 1.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    EndPrimitive();

//  draw m1
    tDir = normalize(m1_gs[1]);
    gl_Position = vec4(pos, 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(1.0, 0.0, 0.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    tDir = normalize(m1_gs[1]);
    gl_Position = vec4(pos + len * tDir, 1.0);
    eDir = normalize(eyePos - vert_gs[1]);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(1.0, 0.0, 0.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    EndPrimitive();

//  draw m2
    tDir = normalize(m2_gs[1]);
    gl_Position = vec4(pos, 1.0);
    eDir = normalize(eyePos - gl_Position.xyz);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 0.0, 1.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    tDir = normalize(m2_gs[1]);
    gl_Position = vec4(pos + len * tDir, 1.0);
    eDir = normalize(eyePos - vert_gs[1]);
    lDir = normalize(light.position - gl_Position).xyz;
    color = vec4(0.0, 0.0, 1.0, 1.0);
    gl_Position = MVP * gl_Position;
    EmitVertex();

    EndPrimitive();
}
