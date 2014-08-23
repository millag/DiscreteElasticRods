#include "Scene.h"
#include "Utils.h"
#include "HairGenerator.h"

RenderObject* createBall(const Mesh* mesh);
Hair* createHair(const RenderObject* object);

//===================================== Scene ===========================================

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
    m_meshes.reserve(10);
    m_renderObjects.reserve(10);

//    initialize scene bounding volume
    mg::Real size = 30;
    m_boundingVolume.reshape(mg::Vec3D(-size, -size, -size), mg::Vec3D(size, size, size));

    unsigned meshId = m_meshes.size();
    m_meshMap["ball"] = meshId;
    Mesh* mesh = Mesh::createSphere(meshId);
    m_meshes.push_back(mesh);

    RenderObject* ball = createBall(mesh);
    m_renderObjects.push_back(ball);

    m_hair = createHair(ball);

    m_spiral = new Spiral();
    m_spiral->init(ball);
}


void Scene::update(mg::Real dt)
{
//    m_hair->update(dt);
    m_spiral->update(dt);
}

void Scene::findObjectsWithinDistance(const mg::Vec3D& pos, mg::Real dist, std::vector<RenderObject*>& o_objects)
{
    typedef std::vector<RenderObject*>::const_iterator Iter;
    for (Iter it = m_renderObjects.begin(); it != m_renderObjects.end(); ++it)
    {
        if ((pos - (*it)->getPosition()).length() < dist + (*it)->getBoundingRadius())
        {
            o_objects.push_back(*it);
        }
    }
}

//====================================== utility functions =============================================

RenderObject *createBall(const Mesh* mesh)
{
    mg::Matrix4D transform;
    mg::matrix_uniform_scale(transform, (mg::Real)1);
    mg::matrix_set_translation(transform, (mg::Real)0, (mg::Real)0, (mg::Real)0);

    return new RenderObject(mesh, transform, -1);
}

Hair* createHair(const RenderObject* object)
{
    std::vector<unsigned> fidx(40);
    for (unsigned i = 0; i < fidx.size(); ++i)
    {
        fidx[i] = 2 * i;
    }

    Hair* hair = new Hair();
    HairGenerator::generateCurlyHair(object, fidx, *hair);

    return hair;
}
