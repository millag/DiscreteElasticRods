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

layout (lines_adjacency) in;
layout (line_strip, max_vertices = 50) out;

in vec3 wvert_cs[];
in vec2 uv_cs[];
in vec3 wnormal_cs[];

out vec3 eDir;
out vec3 lDir;
out vec3 tDir;
out vec2 uv_fr;
out vec3 normal_fr;

vec3 catmull_rom(float t)
{
    return 0.5f * ((2 * wvert_cs[1]) +
                  (wvert_cs[2] - wvert_cs[0]) * t +
                  (2 * wvert_cs[0] - 5 * wvert_cs[1] + 4 * wvert_cs[2] - wvert_cs[3]) * t * t +
                  (3 * wvert_cs[1] - wvert_cs[0] - 3 * wvert_cs[2] + wvert_cs[3]) * t * t * t);
}

vec3 catmull_rom_normal(float t)
{
    return 0.5f * ((2 * wnormal_cs[1]) +
                  (wnormal_cs[2] - wnormal_cs[0]) * t +
                  (2 * wnormal_cs[0] - 5 * wnormal_cs[1] + 4 * wnormal_cs[2] - wnormal_cs[3]) * t * t +
                  (3 * wnormal_cs[1] - wnormal_cs[0] - 3 * wnormal_cs[2] + wnormal_cs[3]) * t * t * t);
}

void main ()
{
    const int nVertices = 2;
    vec3 pos[nVertices];
    vec3 norm[nVertices];

    float step = 1.0 / (pos.length() - 1);
    for (int i = 0; i < pos.length(); i++)
    {
        pos[i] = catmull_rom(i * step);
        norm[i] = wnormal_cs[1];// normalize(catmull_rom_normal(i * step));
    }

    int prev, next;
    for (int i = 0; i < pos.length(); i++)
    {
        gl_Position = vec4(pos[i], 1.0);

        eDir = normalize(eyePos - pos[i]);
        lDir = normalize(light.position - gl_Position).xyz;
        prev = clamp(i - 1, 0, pos.length() - 1);
        next = clamp(i + 1, 0, pos.length() - 1);
        tDir = normalize(pos[next] - pos[prev]);

        gl_Position = MVP * gl_Position;
        EmitVertex();
    }
    EndPrimitive();

    for (int i = 0; i < pos.length(); i++)
    {
        gl_Position = vec4(pos[i], 1.0);

        eDir = normalize(eyePos - pos[i]);
        lDir = normalize(light.position - gl_Position).xyz;
        tDir = norm[i];

        gl_Position = MVP * gl_Position;
        EmitVertex();

        gl_Position = vec4(pos[i], 1.0) + 0.1 * vec4(norm[i], 0.0);

        eDir = normalize(eyePos - pos[i]);
        lDir = normalize(light.position - gl_Position).xyz;
        tDir = norm[i];

        gl_Position = MVP * gl_Position;
        EmitVertex();

        EndPrimitive();
    }
}
