#include "Scene.h"
#include "ngl/Util.h"
#include "Utils.h"
#include "HairyObject.h"

Mesh* createSphere(int id);
HairyObject* createHairyBall(const Mesh* mesh);
RenderObject* createBall(const Mesh* mesh);

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

//    delete hair objects
    typedef std::vector<Hair*>::const_iterator HIter;
    for (HIter it = m_hairs.begin(); it != m_hairs.end(); ++it)
    {
        delete (*it);
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
    m_hairs.reserve(10);

//    initialize scene bounding volume
    ngl::Real size = 30;
    m_boundingVolume.reshape(ngl::Vec4(-size, -size, -size), ngl::Vec4(size, size, size));

    unsigned meshId = m_meshes.size();
    m_meshMap["ball"] = meshId;
    Mesh* mesh = createSphere(meshId);
    m_meshes.push_back(mesh);

    RenderObject* ball = createBall(mesh);
    m_renderObjects.push_back(ball);

//    m_hair = new HairStrand();
//    m_hair->initialize(ball);

    m_spiral = new Spiral();
    m_spiral->init(ball);

//    HairyObject* hairyBall = createHairyBall(mesh);
//    m_renderObjects.push_back(hairyBall);


//    Hair* hair = new Hair(this);
//    m_hairs.push_back(hair);
//    hairyBall->attachHair(hair, 1);
}


void Scene::update(mg::Real dt)
{
//    typedef std::vector<Hair*>::const_iterator Iter;
//    for (Iter it = m_hairs.begin(); it != m_hairs.end(); ++it)
//    {
//        (*it)->update(dt);
//    }

//    m_hair->update(dt);

    m_spiral->update(dt);
}

void Scene::findObjectsWithinDistance(const ngl::Vec4& pos, ngl::Real dist, std::vector<RenderObject*>& o_objects)
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


Mesh* createSphere(int id)
{
    Mesh* mesh = new Mesh(id);

    const unsigned divu = 20;
    const unsigned divv = 10;
    const ngl::Real radius = 1.0;

    mesh->m_vertices.push_back(ngl::Vec4(0,radius,0));
    for (unsigned i = 1; i < divv; i++)
    {
        ngl::Real y = std::cos( ((float)i / divv) * ngl::PI ) * radius;
        ngl::Real r = std::sin( ((float)i / divv) * ngl::PI ) * radius;

        for (unsigned j = 0; j < divu; j++)
        {
            ngl::Real x = std::cos( ((float)j / divu) * ngl::PI * 2) * r;
            ngl::Real z = std::sin( ((float)j / divu) * ngl::PI * 2) * r;
            mesh->m_vertices.push_back(ngl::Vec4(x,y,z));
        }
    }
    mesh->m_vertices.push_back(ngl::Vec4(0,-radius,0));

    for (unsigned j = 0; j < divu; j++)
    {
        unsigned idxp1 = 0;
        unsigned idxc1 = j % divu + 1;
        unsigned idxc2 = (j + 1) % divu + 1;

        mesh->m_vindices.push_back(idxp1);
        mesh->m_vindices.push_back(idxc2);
        mesh->m_vindices.push_back(idxc1);

        idxc1 = (divv - 2) * divu + j % divu + 1;
        idxc2 = (divv - 2) * divu + (j + 1) % divu + 1;
        idxp1 = (divv - 1) * divu + 1;

        mesh->m_vindices.push_back(idxc1);
        mesh->m_vindices.push_back(idxp1);
        mesh->m_vindices.push_back(idxc2);

    }

    for (unsigned i = 1; i < divv - 1; i++)
    {
        for (unsigned j = 0; j < divu; j++)
        {
            unsigned idxp1 = (i - 1) * divu + j % divu + 1;
            unsigned idxp2 = (i - 1) * divu + (j + 1) % divu + 1;
            unsigned idxc1 = i * divu + j % divu + 1;
            unsigned idxc2 = i * divu + (j + 1) % divu + 1;

            mesh->m_vindices.push_back(idxp1);
            mesh->m_vindices.push_back(idxp2);
            mesh->m_vindices.push_back(idxc1);


            mesh->m_vindices.push_back(idxc2);
            mesh->m_vindices.push_back(idxc1);
            mesh->m_vindices.push_back(idxp2);
        }
    }

    Mesh::calcNormals( *mesh,  Mesh::FLAT );
    return mesh;
}

HairyObject *createHairyBall(const Mesh* mesh)
{
    ngl::Real s = 0.85;
    ngl::Vec4 tr = ngl::Vec4(0, 1.5, 0);
    ngl::Mat4 t;
    t.scale(s, s, s);
    t.translate(tr.m_x, tr.m_y, tr.m_z);

    return new HairyObject(mesh, t, -1);
}

RenderObject *createBall(const Mesh* mesh)
{
    ngl::Real s = 0.1;
    ngl::Vec4 tr = ngl::Vec4(-1.5, 2, 0);
    ngl::Mat4 t;
    t.scale(s, s, s);
    t.translate(tr.m_x, tr.m_y, tr.m_z);

    return new RenderObject(mesh, t, -1);
}
