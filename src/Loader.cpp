#include "Loader.h"

#include "ObjLoader.h"
#include "HairGenerator.h"

struct Loader::PImpl
{
    Mesh* loadMesh(const char *fileName);
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

    Mesh* mesh = m_impl->loadMesh("assets/medusa1.obj");
    unsigned meshId = scene->m_meshes.size();
    mesh->setId(meshId);
    scene->m_meshes.push_back(mesh);

    mg::Matrix4D transform;
    transform.identity();
    RenderObject* object = new RenderObject(mesh, transform);
    scene->m_renderObjects.push_back(object);

    mg::Matrix4D shape, rot;
    mg::matrix_rotation_euler(shape, mg::rad((mg::Real)-23.853), (mg::Real)0, (mg::Real)0, mg::euler_order_xyz);
    mg::matrix_scale(rot, (mg::Real)1.967, (mg::Real)2.953, (mg::Real)2.355);
    shape = shape * rot;
    mg::matrix_set_translation(shape, (mg::Real)0, (mg::Real)2.664, (mg::Real)0);
    object->addCollisionShape(shape);

    unsigned nFidx = 224;
    unsigned fidx[] = {861, 862, 863, 864, 875, 876, 877, 878, 879, 880, 881, 882, 883, 884, 891, 892, 893, 894, 895, 896,
                       897, 898, 899, 900, 901, 902, 903, 904, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920,
                       1336, 1337, 1338, 1339, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355, 1356,
                       1357, 1358, 1359, 1363, 1364, 1365, 1366, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375, 1376,
                       1377, 1378, 1379, 1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1391, 1392, 1393,
                       1394, 1395, 1396, 1397, 1398, 1399, 1400, 1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1413,
                       1416, 1417, 1420, 1421, 3334, 3335, 3336, 3337, 3348, 3349, 3350, 3351, 3352, 3353, 3354, 3355, 3356,
                       3357, 3364, 3365, 3366, 3367, 3368, 3369, 3370, 3371, 3372, 3373, 3374, 3375, 3376, 3377, 3382, 3383,
                       3384, 3385, 3386, 3387, 3388, 3389, 3390, 3391, 3392, 3393, 3804, 3805, 3806, 3807, 3812, 3813, 3814,
                       3815, 3816, 3817, 3818, 3819, 3820, 3821, 3822, 3823, 3824, 3825, 3826, 3827, 3831, 3832, 3833, 3834,
                       3835, 3836, 3837, 3838, 3839, 3840, 3841, 3842, 3843, 3844, 3845, 3846, 3847, 3848, 3849, 3850, 3851,
                       3852, 3853, 3854, 3855, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3863, 3864, 3865, 3866, 3867, 3868,
                       3869, 3870, 3872, 3873, 3874, 3875, 3876, 3877, 3878, 3879, 3880, 3881, 3882, 3883};
    std::vector<unsigned> findices(nFidx);
    std::copy(fidx, fidx + nFidx, findices.begin());

    scene->m_hair = new Hair();
    HairGenerator::generateCurlyHair(object, findices, *scene->m_hair);

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


Mesh *Loader::PImpl::loadMesh(const char* fileName)
{
    ObjLoader loader;
    if (!loader.loadFile(fileName))
    {
        std::cerr << "Error loading mesh from file " << fileName << std::endl;
        return NULL;
    }

    Mesh* mesh = new Mesh();
    mesh->m_vertices.resize(loader.getVertices().size());
    std::copy(loader.getVertices().begin(), loader.getVertices().end(), mesh->m_vertices.begin());
    mesh->m_vindices.resize(loader.getVIndices().size());
    std::copy(loader.getVIndices().begin(), loader.getVIndices().end(), mesh->m_vindices.begin());

    if (loader.getNormals().size() != loader.getVertices().size() || loader.getVNIndices().size() != loader.getVIndices().size())
    {
        Mesh::computeNormals(*mesh, Mesh::ShadingMode::GOURAUD);
    } else
    {
        mesh->m_normals.resize(loader.getNormals().size());
        std::copy(loader.getNormals().begin(), loader.getNormals().end(), mesh->m_normals.begin());
    }

    return mesh;
}


