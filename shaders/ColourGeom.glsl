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

uniform mat4 M;
uniform mat4 MV;
uniform mat4 MVP;
// array of lights
uniform Lights light;
uniform vec3 eyePos;

layout (lines_adjacency) in;
layout (line_strip, max_vertices = 10) out;

out vec4 eDir;
out vec4 lDir;
out vec4 tDir;

vec4 catmull_rom(float t)
{
    return 0.5f * ((2 * gl_in[1].gl_Position) +
                  (gl_in[2].gl_Position - gl_in[0].gl_Position) * t +
                  (2 * gl_in[0].gl_Position - 5 * gl_in[1].gl_Position + 4 * gl_in[2].gl_Position - gl_in[3].gl_Position) * t * t +
                  (3 * gl_in[1].gl_Position - gl_in[0].gl_Position - 3 * gl_in[2].gl_Position + gl_in[3].gl_Position) * t * t * t);
}

void main ()
{
    const int nVertices = 6;
    vec4 pos[nVertices];

    float step = 1.0 / (pos.length() - 1);
    for (int i = 0; i < pos.length(); i++)
    {
        pos[i] = catmull_rom(i * step);
    }

    int prev, next;
    vec4 viewPos, worldPos;
    for (int i = 0; i < pos.length(); i++)
    {
        gl_Position = pos[i];

        worldPos = M * gl_Position;
        viewPos = MV * gl_Position;
        eDir = normalize(vec4(eyePos, 1.0) - worldPos);
        lDir = normalize(light.position - worldPos);
        prev = clamp(i - 1, 0, pos.length() - 1);
        next = clamp(i + 1, 0, pos.length() - 1);
        tDir =  M * normalize(pos[next] - pos[prev]);

        gl_Position = MVP * gl_Position;
        EmitVertex();
    }
    EndPrimitive();
}
