#include "Loader.h"

#include <iostream>
#include <fstream>
#include "HairGenerator.h"

struct Loader::PImpl
{
    void loadMesh(const char *fileName, Mesh &o_mesh);
};

Loader::Loader()
{
    m_impl = new PImpl();
}

Loader::~Loader()
{
    delete m_impl;
}

Scene* Loader::loadScene(const char*)
{
    Scene* scene = new Scene();

    mg::Real size = 30;
    scene->m_boundingVolume.reshape(mg::Vec3D(-size, -size, -size), mg::Vec3D(size, size, size));

//    unsigned meshId = scene->m_meshes.size();
//    Mesh* mesh = m_impl->load_obj("assets/shape1.obj", meshId);
////    m_impl->loadMesh("assets/shape1.obj", *mesh);
////    Mesh* mesh = new Mesh(meshId, Mesh::PrimitiveMode::TRIANGLES);
////    m_impl->loadMesh("assets/shape1.obj", *mesh);
//    scene->m_meshes.push_back(mesh);

//    mg::Matrix4D transform;
//    transform.identity();
//    RenderObject* object = new RenderObject(mesh, transform);
//    scene->m_renderObjects.push_back(object);

//    mg::Matrix4D shape, rot;
//    mg::matrix_rotation_euler(shape, mg::rad((mg::Real)-23.853), (mg::Real)0, (mg::Real)0, mg::euler_order_xyz);
//    mg::matrix_scale(rot, (mg::Real)1.967, (mg::Real)2.953, (mg::Real)2.355);
//    shape = shape * rot;
//    mg::matrix_set_translation(shape, (mg::Real)0, (mg::Real)2.664, (mg::Real)0);
//    object->addCollisionShape(shape);

//    unsigned nFidx = 4;
//    unsigned fidx[] = {0, 1, 2, 3};
//    std::vector<unsigned> findices(nFidx);
//    std::copy(fidx, fidx + nFidx, findices.begin());

//    scene->m_hair = new Hair();
//    HairGenerator::generateCurlyHair(object, findices, *scene->m_hair);

    return scene;
}


Scene* Loader::loadTestScene()
{
    Scene* scene = new Scene();

    mg::Real size = 30;
    scene->m_boundingVolume.reshape(mg::Vec3D(-size, -size, -size), mg::Vec3D(size, size, size));

//    create mesh
    unsigned meshId = scene->m_meshes.size();
    Mesh* mesh = Mesh::createSphere(meshId);
    scene->m_meshes.push_back(mesh);
//    create render object
    mg::Matrix4D transform;
    mg::matrix_uniform_scale(transform, (mg::Real)1);
    mg::matrix_set_translation(transform, (mg::Real)0, (mg::Real)0, (mg::Real)0);

    RenderObject* object = new RenderObject(mesh, transform, -1);
    scene->m_renderObjects.push_back(object);

    std::vector<unsigned> fidx(120);
    for (unsigned i = 0; i < fidx.size(); ++i)
    {
        fidx[i] = 2 * i;
    }
//    create hair object
    scene->m_hair = new Hair();
//    generate hair
    HairGenerator::generateCurlyHair(object, fidx, *scene->m_hair);

//    add collision shape for the object
    mg::Real radius = object->getMeshBoundingRadius() + scene->m_hair->m_params->m_thickness;
    mg::Matrix4D shape;
    mg::matrix_scale(shape, radius, radius, radius);
    mg::matrix_set_translation(shape, object->getMeshAABB().getCenter());
    object->addCollisionShape(shape);

    scene->m_spiral = new Spiral();
    scene->m_spiral->init(object);

    return scene;
}


void Loader::PImpl::loadMesh(const char* fileName, Mesh& o_mesh)
{

//    Mesh::computeNormals(o_mesh, Mesh::ShadingMode::GOURAUD);
}


