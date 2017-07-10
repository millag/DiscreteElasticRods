#version 410

/// @brief the vertex passed in
layout (location = 0) in vec3 inVert;
/// @brief the normal passed in
layout (location = 2) in vec3 inNormal;
/// @brief the in uv
layout (location = 1) in vec2 inUV;
/// @brief flag to indicate if model has unit normals if not normalize
uniform bool Normalize;
// the eye position of the camera
uniform vec3 viewerPos;
/// @brief the current fragment normal for the vert being processed
out vec3 fragmentNormal;


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
// our material
uniform Materials material;
// array of lights
uniform Lights light;
// direction of the lights used for shading
out vec3 lightDir;
// out the blinn half vector
out vec3 halfVector;
out vec3 eyeDirection;
out vec3 vPosition;

uniform mat4 MV;
uniform mat4 MVP;
uniform mat3 normalMatrix;
uniform mat4 M;


void main()
{
	// calculate the fragments surface normal
	fragmentNormal = (normalMatrix*inNormal);


	if (Normalize == true)
	{
		fragmentNormal = normalize(fragmentNormal);
	}
	// calculate the vertex position
	gl_Position = MVP*vec4(inVert,1.0);

	vec4 worldPosition = M * vec4(inVert, 1.0);
	eyeDirection = normalize(viewerPos - worldPosition.xyz);
	// Get vertex position in eye coordinates
	// Transform the vertex to eye co-ordinates for frag shader
	/// @brief the vertex in eye co-ordinates  homogeneous
	vec4 eyeCord=MV*vec4(inVert,1);

	vPosition = eyeCord.xyz / eyeCord.w;;

	float dist;

	lightDir=vec3(light.position.xyz-eyeCord.xyz);
	dist = length(lightDir);
	lightDir/= dist;
	halfVector = normalize(eyeDirection + lightDir);

}
