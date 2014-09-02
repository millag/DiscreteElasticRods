#include "SceneLoader.h"
#include <map>

#include "ObjLoader.h"
#include "HairGenerator.h"
#include "BasicParser.h"

struct mesh_object
{
    mesh_object():m_id(-1) { }
    long int m_id;
    std::string m_filename;
};

struct object_3D
{
    object_3D():m_id(-1), m_meshId(-1) { m_transform.identity(); }
    long int m_id;
    long int m_meshId;
    mg::Matrix4D m_transform;
    std::vector<mg::Matrix4D> m_collisionShapes;
};

struct hair_object
{
    hair_object():m_id(-1), m_objId(-1), m_faceCnt(0) { }
    long int m_id;
    long int m_objId;
    std::string m_type;
    unsigned m_faceCnt;
    std::vector<unsigned> m_faceList;
};

struct hair_params
{
    hair_params():m_id(-1), m_hairId(-1) { }
    long int m_id;
    long int m_hairId;
//    TODO: add all other props
};

struct scene_object
{
    std::map<unsigned, mesh_object> m_meshMap;
    std::map<unsigned, object_3D> m_object3DMap;
    std::map<unsigned, hair_object> m_hairMap;
    std::map<unsigned, hair_params> m_hairParamsMap;
};




struct SceneLoader::PImpl
{
public:
    inline std::string &stripComment(std::string& s);
    bool parseScene(std::ifstream& ifs, scene_object& o_scene);
    bool parseMesh(std::ifstream& ifs, scene_object& o_scene);
    bool parseObject3D(std::ifstream& ifs, scene_object& o_scene);
    bool parseHairObject(std::ifstream& ifs, scene_object& o_scene);
    bool parseHairParams(std::ifstream& ifs, scene_object& o_scene);

    unsigned m_lineNumber;
};




SceneLoader::SceneLoader()
{
    m_impl = new PImpl();
}

SceneLoader::~SceneLoader()
{
    delete m_impl;
}

Scene* SceneLoader::loadTestScene()
{
    Scene* scene = new Scene();

    mg::Real size = 30;
    scene->m_boundingVolume.reshape(mg::Vec3D(-size, -size, -size), mg::Vec3D(size, size, size));

//    create mesh
    unsigned meshId = scene->m_meshes.size();
    Mesh* mesh = Mesh::createSphere(meshId, 20, 10);
    scene->m_meshes.push_back(mesh);
//    create render object
    mg::Matrix4D transform;
    mg::matrix_uniform_scale(transform, (mg::Real)1);
    mg::matrix_set_translation(transform, (mg::Real)0, (mg::Real)0, (mg::Real)0);

    RenderObject* object = new RenderObject(mesh, transform, -1);
    scene->m_renderObjects.push_back(object);

    std::vector<unsigned> fidx(object->getMesh()->getNPrimitives() / 2);
    for (unsigned i = 0; i < fidx.size(); ++i)
    {
        fidx[i] = 2 * i;
    }
//    create hair object
    Hair* hair = new Hair();
    scene->m_hairs.push_back(hair);
//    generate hair
    HairGenerator::generateCurlyHair(object, fidx, *hair);

//    add collision shape for the object
    mg::Real radius = object->getMeshBoundingRadius() + hair->m_params->m_thickness;
    mg::Matrix4D shape;
    mg::matrix_scale(shape, radius, radius, radius);
    mg::matrix_set_translation(shape, object->getMeshAABB().getCenter());
    object->addCollisionShape(shape);

    scene->m_spiral = new Spiral();
    scene->m_spiral->init(object);

    return scene;
}


Scene* SceneLoader::loadScene(const char* filename)
{

    if (filename == NULL)
        return false;

    std::ifstream ifs;
    ifs.open(filename, std::ios::in);
    if (ifs.fail())
    {
        std::cerr << "Can't open file " << filename << std::endl;
        ifs.close();
        return NULL;
    }

    scene_object scene_obj;
    bool success = m_impl->parseScene(ifs, scene_obj);
    if (!success)
    {
        std::cerr << "Error parsing file " << filename  <<
                     "\nError at line: " << m_impl->m_lineNumber << std::endl;
        ifs.close();
        return NULL;
    }

//    construct scene
    Scene* scene = new Scene();

    mg::Real size = 30;
    scene->m_boundingVolume.reshape(mg::Vec3D(-size, -size, -size), mg::Vec3D(size, size, size));

//    create meshes
    typedef std::map<unsigned, mesh_object>::iterator MIter;
    for (MIter it = scene_obj.m_meshMap.begin(); it != scene_obj.m_meshMap.end(); ++it)
    {
        Mesh* mesh = loadMesh(it->second.m_filename.c_str());
        if (mesh == NULL)
        {
            it->second.m_id = -1;
            continue;
        }
        it->second.m_id = scene->m_meshes.size();
        mesh->setId(it->second.m_id);
        scene->m_meshes.push_back(mesh);
    }
//    create render objects
    typedef std::map<unsigned, object_3D>::iterator OIter;
    for (OIter it = scene_obj.m_object3DMap.begin(); it != scene_obj.m_object3DMap.end(); ++it)
    {
        if (!scene_obj.m_meshMap.count(it->second.m_meshId) || scene_obj.m_meshMap[it->second.m_meshId].m_id < 0)
        {
            it->second.m_id = -1;
            continue;
        }
        unsigned idx = scene_obj.m_meshMap[it->second.m_meshId].m_id;
        RenderObject* object = new RenderObject(scene->m_meshes[idx], it->second.m_transform);
        for (unsigned i = 0; i < it->second.m_collisionShapes.size(); ++i)
        {
            object->addCollisionShape(it->second.m_collisionShapes[i]);
        }
        it->second.m_id = scene->m_renderObjects.size();
        scene->m_renderObjects.push_back(object);
    }
//    create hair objects
    typedef std::map<unsigned, hair_object>::iterator Iter;
    for (Iter it = scene_obj.m_hairMap.begin(); it != scene_obj.m_hairMap.end(); ++it)
    {
        if (!scene_obj.m_object3DMap.count(it->second.m_objId) || scene_obj.m_object3DMap[it->second.m_objId].m_id < 0)
        {
            it->second.m_id = -1;
            continue;
        }

        unsigned idx = scene_obj.m_object3DMap[it->second.m_objId].m_id;
        Hair* hair = new Hair();
        it->second.m_id = scene->m_hairs.size();
        scene->m_hairs.push_back(hair);
        if (it->second.m_type.compare("straight"))
        {
            HairGenerator::generateCurlyHair(scene->m_renderObjects[idx], it->second.m_faceList, *hair);
        } else
        {
            HairGenerator::generateStraightHair(scene->m_renderObjects[idx], it->second.m_faceList, *hair);
        }
    }

    ifs.close();
    return scene;
}


Mesh *SceneLoader::loadMesh(const char* fileName)
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


inline std::string &SceneLoader::PImpl::stripComment(std::string& s)
{
    std::size_t found = s.find('#');
    if (found == std::string::npos)
        return s;
    return s.erase(found, s.size());
}

bool SceneLoader::PImpl::parseScene(std::ifstream &ifs, scene_object &o_scene)
{
    bool success = true;

    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, token, '\n');
        ++m_lineNumber;

        stripComment(token);
        BasicParser::trim(token);

        if (token.empty())
        {
            continue;
        }

        if (token.compare("mesh") == 0)
        {
            success = parseMesh(ifs, o_scene);
            if (!success)
            {
                break;
            }
            continue;
        }

        if (token.compare("object3D") == 0)
        {
            success = parseObject3D(ifs, o_scene);
            if (!success)
            {
                break;
            }
            continue;
        }

        if (token.compare("hair") == 0)
        {
            success = parseHairObject(ifs, o_scene);
            if (!success)
            {
                break;
            }
            continue;
        }

        if (token.compare("hairParams") == 0)
        {
            success = parseHairParams(ifs, o_scene);
            if (!success)
            {
                break;
            }
            continue;
        }

        success = false;
        std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
        break;
    }

    return success;
}

bool SceneLoader::PImpl::parseMesh(std::ifstream &ifs, scene_object &o_scene)
{
    bool success = false;

    mesh_object mesh;

    std::string line;
    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line, '\n');
        ++m_lineNumber;

        stripComment(line);
        BasicParser::trim(line);
        if (line.empty())
        {
            continue;
        }

        token = BasicParser::parseToken(line);
        if (token.compare("id") == 0)
        {
            line.erase(0, token.length());
            mesh.m_id = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("filename") == 0)
        {
            line.erase(0, token.length());
            mesh.m_filename = BasicParser::trim(line);
            continue;
        }

        if (token.compare("/mesh") == 0)
        {
            success = true;
            break;
        }

        success = false;
        std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
        break;

    }

    if (!success || mesh.m_id < 0 || o_scene.m_meshMap.count(mesh.m_id))
    {
        success = false;
    } else
    {
        o_scene.m_meshMap[mesh.m_id] = mesh;
    }
    return success;
}


bool SceneLoader::PImpl::parseObject3D(std::ifstream &ifs, scene_object &o_scene)
{
    bool success = false;

    object_3D object;

    std::string line;
    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line, '\n');
        ++m_lineNumber;

        stripComment(line);
        BasicParser::trim(line);
        if (line.empty())
        {
            continue;
        }

        token = BasicParser::parseToken(line);
        if (token.compare("id") == 0)
        {
            line.erase(0, token.length());
            object.m_id = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("meshId") == 0)
        {
            line.erase(0, token.length());
            object.m_meshId = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("transform") == 0)
        {
            line.erase(0, token.length());
            BasicParser::parseMatrix4D(line, object.m_transform);
            continue;
        }

        if (token.compare("collisionShape") == 0)
        {
            line.erase(0, token.length());
            unsigned idx = object.m_collisionShapes.size();
            object.m_collisionShapes.resize(idx + 1);
            BasicParser::parseMatrix4D(line, object.m_collisionShapes[idx]);
            continue;
        }

        if (token.compare("/object3D") == 0)
        {
            success = true;
            break;
        }

        success = false;
        std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
        break;

    }

    if (!success || object.m_id < 0 || object.m_meshId < 0 || o_scene.m_object3DMap.count(object.m_id))
    {
        success = false;
    } else
    {
        o_scene.m_object3DMap[object.m_id] = object;
    }
    return success;
}

bool SceneLoader::PImpl::parseHairObject(std::ifstream& ifs, scene_object& o_scene)
{
    bool success = false;

    hair_object object;

    std::string line;
    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line, '\n');
        ++m_lineNumber;

        stripComment(line);
        BasicParser::trim(line);
        if (line.empty())
        {
            continue;
        }

        token = BasicParser::parseToken(line);
        if (token.compare("id") == 0)
        {
            line.erase(0, token.length());
            object.m_id = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("objId") == 0)
        {
            line.erase(0, token.length());
            object.m_objId = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("type") == 0)
        {
            line.erase(0, token.length());
            BasicParser::ltrim(line);
            object.m_type = BasicParser::parseWord(line);
            continue;
        }

        if (token.compare("faceCnt") == 0)
        {
            line.erase(0, token.length());
            object.m_faceCnt = BasicParser::parseInt(line);
            object.m_faceList.resize(object.m_faceCnt);
            continue;
        }

        if (token.compare("faceList") == 0)
        {
            line.erase(0, token.length());
            BasicParser::parseIntList(line, object.m_faceList);
            continue;
        }

        if (token.compare("/hair") == 0)
        {
            success = true;
            break;
        }

        success = false;
        std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
        break;

    }

    if (!success || object.m_id < 0 || object.m_objId < 0 || o_scene.m_hairMap.count(object.m_id))
    {
        success = false;
    } else
    {
        o_scene.m_hairMap[object.m_id] = object;
    }
    return success;
}

bool SceneLoader::PImpl::parseHairParams(std::ifstream& ifs, scene_object& o_scene)
{
    bool success = false;

    hair_params object;

    std::string line;
    std::string token;
    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line, '\n');
        ++m_lineNumber;

        stripComment(line);
        BasicParser::trim(line);
        if (line.empty())
        {
            continue;
        }

        token = BasicParser::parseToken(line);
        if (token.compare("id") == 0)
        {
            line.erase(0, token.length());
            object.m_id = BasicParser::parseInt(line);
            continue;
        }

        if (token.compare("hairId") == 0)
        {
            line.erase(0, token.length());
            object.m_hairId = BasicParser::parseInt(line);
            continue;
        }

//        if (token.compare("type") == 0)
//        {
//            line.erase(0, token.length());
//            object.m_type = BasicParser::parseWord(line);
//            continue;
//        }

//        if (token.compare("faceCnt") == 0)
//        {
//            line.erase(0, token.length());
//            object.m_faceCnt = BasicParser::parseInt(line);
//            continue;
//        }

//        if (token.compare("faceList") == 0)
//        {
//            line.erase(0, token.length());
//            object.m_faceCnt = BasicParser::parseInt(line);
//            continue;
//        }

        if (token.compare("/hairParams") == 0)
        {
            success = true;
            break;
        }

        success = false;
        std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
        break;

    }

    if (!success || object.m_id < 0 || object.m_hairId < 0 || o_scene.m_hairParamsMap.count(object.m_id))
    {
        success = false;
    } else
    {
        o_scene.m_hairParamsMap[object.m_id] = object;
    }
    return success;
}
