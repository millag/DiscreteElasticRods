#version 410

uniform vec3 color;

/// output fragment colour
layout (location = 0) out vec4 fragColor;

void main(void)
{
	fragColor = vec4( color, 1.f );
}
