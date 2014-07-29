#version 410
/// @file Colour.fs
/// @brief a basic unshaded solid colour shader
/// @brief the colour to shade draw with
uniform vec4 Colour;
in vec4 eDir;
in vec4 lDir;
in vec4 tDir;

out vec4 fragColour;

void main()
{
    vec4 e = normalize(eDir);
    vec4 l = normalize(lDir);
    vec4 t = normalize(tDir);

    float dotTL = dot(l, t);
    float dotTE = dot(t, e);
    vec4 diffuse = Colour * sqrt(1.0 - dotTL * dotTL);
    vec4 specular = vec4(1.,1.,1.,1.) *  pow((dotTL * dotTE + sqrt(1.0 - dotTL * dotTL) * sqrt(1.0 - dotTE * dotTE)), 0.2);
    fragColour = diffuse;
}

