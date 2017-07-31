#version 410

uniform vec4 color;

/// output fragment colour
layout (location = 0) out vec4 fragColor;

void main(void)
{
	fragColor = color;
}
