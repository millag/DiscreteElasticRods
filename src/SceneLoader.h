#ifndef SCENELOADER_H
#define SCENELOADER_H

#include "Scene.h"

#include <memory>


class SceneLoader
{
public:
	SceneLoader();
	~SceneLoader();

	bool loadTestScene( Scene& scene );
	bool loadScene( const char *filename, Scene& scene );
	bool loadMesh( const char *fileName, Mesh& mesh );

private:
	struct PImpl;
	std::unique_ptr<PImpl> m_impl;
};

#endif // SCENELOADER_H
