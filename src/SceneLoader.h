#pragma once

#include "Scene.h"

class SceneLoader
{
public:
	SceneLoader();
	~SceneLoader();

	bool loadTestScene( Scene& scene );
	bool loadScene( const char *filename, Scene& scene );
	bool loadMesh( const char *fileName, Mesh& mesh );

private:
	NON_COPYABLE( SceneLoader );

	struct PImpl;
	std::unique_ptr<PImpl> m_impl;
};
