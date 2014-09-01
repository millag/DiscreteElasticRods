#include "Scene.h"

Scene::Scene():m_spiral(NULL) { }

Scene::~Scene()
{
    if (m_spiral != NULL)
    {
        delete m_spiral;
    }

//    delete hair objects
    typedef std::vector<Hair*>::const_iterator Iter;
    for (Iter it = m_hairs.begin(); it != m_hairs.end(); ++it)
    {
        delete (*it);
    }

//    delete render objects
    typedef std::vector<RenderObject*>::const_iterator RIter;
    for (RIter it = m_renderObjects.begin(); it != m_renderObjects.end(); ++it)
    {
        delete (*it);
    }

//    delete meshes
    typedef std::vector<Mesh*>::const_iterator MIter;
    for (MIter it = m_meshes.begin(); it != m_meshes.end(); ++it)
    {
        delete (*it);
    }
}

void Scene::initialize()
{
//    TODO:implement
}


void Scene::update(mg::Real dt)
{
    m_hairs[0]->update(dt);
//    m_spiral->update(dt);
}
