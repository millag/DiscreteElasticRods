#ifndef LOADER_H
#define LOADER_H

#include "Scene.h"

class Loader
{
public:
    Loader();
    ~Loader();

    Scene* loadScene(const std::string& fileName);

private:

    struct PImpl;
    PImpl* m_impl;
};

#endif // LOADER_H
