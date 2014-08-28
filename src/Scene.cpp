#include "Scene.h"

Scene::Scene(): m_hair(NULL), m_spiral(NULL) { }

Scene::~Scene()
{

//    delete hair
    if (m_hair != NULL)
    {
        delete m_hair;
    }

    if (m_spiral != NULL)
    {
        delete m_spiral;
    }

//    delete render objects
    typedef std::vector<RenderObject*>::const_iterator ROIter;
    for (ROIter it = m_renderObjects.begin(); it != m_renderObjects.end(); ++it)
    {
        delete (*it);
    }

//    delete meshes
    typedef std::vector<Mesh*>::const_iterator GIter;
    for (GIter it = m_meshes.begin(); it != m_meshes.end(); ++it)
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
    m_hair->update(dt);
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
