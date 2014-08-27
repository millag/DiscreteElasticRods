#ifndef LOADER_H
#define LOADER_H

#include "Scene.h"

class Loader
{
public:
    Loader();
    ~Loader();

    Scene* loadScene(const char *);
    Scene* loadTestScene();
private:

    struct PImpl;
    PImpl* m_impl;
};

#endif // LOADER_H
