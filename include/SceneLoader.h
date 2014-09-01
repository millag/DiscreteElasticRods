#ifndef SCENELOADER_H
#define SCENELOADER_H

#include "Scene.h"

class SceneLoader
{
public:

    SceneLoader();
    ~SceneLoader();

    Scene* loadScene(const char *filename);
    Scene* loadTestScene();
    Mesh* loadMesh(const char *fileName);
private:

    struct PImpl;
    PImpl* m_impl;
};

#endif // SCENELOADER_H
