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

void Scene::findObjectsWithinDistance(const mg::Vec3D& pos, mg::Real dist, std::vector<RenderObject*>& o_objects)
{
    typedef std::vector<RenderObject*>::const_iterator Iter;
    for (Iter it = m_renderObjects.begin(); it != m_renderObjects.end(); ++it)
    {
        if ((pos - (*it)->getCenter()).length() < dist + (*it)->getBoundingRadius())
        {
            o_objects.push_back(*it);
        }
    }
}
