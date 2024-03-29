#include "SceneLoader.h"
#include "ObjLoader.h"
#include "HairGenerator.h"
#include "BasicParser.h"
#include <map>

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
	unsigned m_faceCnt;
	std::vector<unsigned> m_faceList;

	double m_length;
	double m_lengthVariance;
	double m_helicalRadius;
	double m_helicalPitch;
	double m_density;
	double m_thickness;
	long int m_nParticles;

	mg::Vec3D m_netForce;
	double m_drag;

	long int m_resolveCollisions;
	long int m_resolveSelfInterations;
	double m_selfInterationDist;
	double m_selfStiction;
	double m_selfRepulsion;

	long int m_pbdIter;
	double m_bendStiffness;
	double m_twistStiffness;
	double m_maxElasticForce;

	ElasticRodParams::MinimizationStrategy m_minimizationStrategy;
	double m_minimizationTolerance;
	long int m_minimizationMaxIter;
};


struct scene_object
{
	std::map<unsigned, mesh_object> m_meshMap;
	std::map<unsigned, object_3D> m_object3DMap;
	std::map<unsigned, hair_object> m_hairMap;
};




struct SceneLoader::PImpl
{
public:
	inline std::string &stripComment(std::string& s);
	bool parseScene(std::ifstream& ifs, scene_object& o_scene);
	bool parseMesh(std::ifstream& ifs, scene_object& o_scene);
	bool parseObject3D(std::ifstream& ifs, scene_object& o_scene);
	bool parseHairObject(std::ifstream& ifs, scene_object& o_scene);

	unsigned m_lineNumber;
};

SceneLoader::SceneLoader():
    m_impl( std::make_unique<PImpl>() )
{ }

SceneLoader::~SceneLoader() = default;

bool SceneLoader::loadTestScene( Scene& scene )
{
	const auto size = 30.f;
	scene.m_boundingVolume.reshape( mg::Vec3D( -size, -size, -size ), mg::Vec3D( size, size, size ) );

//	mesh
	const auto meshId = scene.m_meshes.size();
	auto mesh = Mesh::createSphere(meshId, 20, 10);
	scene.m_meshes.push_back( mesh );
//	render object
	mg::Matrix4D transform;
	mg::matrix_uniform_scale( transform, 1.f );
	mg::matrix_set_translation( transform, 0.f, 0.f, 0.f );

	auto object = std::make_shared<RenderObject>( mesh.get(), transform );
	scene.m_renderObjects.push_back( object );

	std::vector<unsigned> fidx( 60 );
	for ( auto i = 0u; i < fidx.size(); ++i )
	{
		fidx[i] = 2*i;
	}

//	hair object
	auto hair = std::make_shared<Hair>();

	HairGenerator generator;
	generator.generateHair( *object, fidx, *hair );

	scene.m_hairs.push_back( hair );

//	add collision shape for the object
	const auto radius = object->getMeshBoundingRadius() + hair->m_params.m_thickness;
	mg::Matrix4D shape;
	mg::matrix_scale( shape, radius, radius, radius );
	mg::matrix_set_translation( shape, object->getMeshAABB().getCenter() );
	object->addCollisionShape( shape );

	scene.m_spiral = std::make_unique<Spiral>();
	scene.m_spiral->initialize( *object );

	return true;
}

bool SceneLoader::loadScene( const char* filename, Scene& scene )
{
	if ( !filename || !*filename )
	{
		return false;
	}

	std::ifstream ifs;
	ifs.open( filename, std::ios::in );
	if ( ifs.fail() )
	{
		ifs.close();
		std::cerr << "Can't open file " << filename << std::endl;
		return false;
	}

	scene_object scene_obj;
	if ( !m_impl->parseScene( ifs, scene_obj ) )
	{
		ifs.close();
		std::cerr << "Error parsing file " << filename
		          << "\nError at line: " << m_impl->m_lineNumber << std::endl;
		return false;
	}

	const auto size = 30.f;
	scene.m_boundingVolume.reshape( mg::Vec3D( -size, -size, -size ), mg::Vec3D( size, size, size ) );

//	meshes
	for ( auto it = scene_obj.m_meshMap.begin(); it != scene_obj.m_meshMap.end(); ++it )
	{
		auto mesh = std::make_shared<Mesh>();
		if( !loadMesh( it->second.m_filename.c_str(), *mesh ) )
		{
			it->second.m_id = -1;
			continue;
		}

		it->second.m_id = scene.m_meshes.size();
		scene.m_meshes.push_back( mesh );
	}

//	render objects
	for ( auto it = scene_obj.m_object3DMap.begin(); it != scene_obj.m_object3DMap.end(); ++it )
	{
		if ( !scene_obj.m_meshMap.count( it->second.m_meshId ) || scene_obj.m_meshMap[it->second.m_meshId].m_id < 0 )
		{
			it->second.m_id = -1;
			continue;
		}

		const auto idx = scene_obj.m_meshMap[it->second.m_meshId].m_id;
		auto object = std::make_shared<RenderObject>( scene.m_meshes[idx].get(), it->second.m_transform );
		for ( auto i = 0u; i < it->second.m_collisionShapes.size(); ++i )
		{
			object->addCollisionShape( it->second.m_collisionShapes[i] );
		}

		it->second.m_id = scene.m_renderObjects.size();
		scene.m_renderObjects.push_back( object );
	}

//	hair objects
	for ( auto it = scene_obj.m_hairMap.begin(); it != scene_obj.m_hairMap.end(); ++it )
	{
		if ( !scene_obj.m_object3DMap.count( it->second.m_objId ) || scene_obj.m_object3DMap[it->second.m_objId].m_id < 0 )
		{
			it->second.m_id = -1;
			continue;
		}

		const auto idx = scene_obj.m_object3DMap[it->second.m_objId].m_id;
		auto hair = std::make_shared<Hair>();

//		assign hair props
		hair->m_params.m_length = it->second.m_length;
		hair->m_params.m_lengthVariance = it->second.m_lengthVariance;
		hair->m_params.m_helicalRadius = it->second.m_helicalRadius;
		hair->m_params.m_helicalPitch = it->second.m_helicalPitch;
		hair->m_params.m_density = it->second.m_density;
		hair->m_params.m_thickness = it->second.m_thickness;
		hair->m_params.m_nParticles = it->second.m_nParticles;

		hair->m_params.m_gravity = it->second.m_netForce;
		hair->m_params.m_drag = it->second.m_drag;

		hair->m_params.m_resolveCollisions = it->second.m_resolveCollisions > 0;
		hair->m_params.m_resolveSelfInterations = it->second.m_resolveSelfInterations > 0;
		hair->m_params.m_selfInterationDist = it->second.m_selfInterationDist;
		hair->m_params.m_selfStiction = it->second.m_selfStiction;
		hair->m_params.m_selfRepulsion = it->second.m_selfRepulsion;

		hair->m_params.m_pbdIter = it->second.m_pbdIter;
		hair->m_params.m_rodParams.setBendStiffness(it->second.m_bendStiffness);
		hair->m_params.m_rodParams.setTwistStiffness(it->second.m_twistStiffness);
		hair->m_params.m_rodParams.m_maxElasticForce = it->second.m_maxElasticForce;

		hair->m_params.m_rodParams.m_strategy = it->second.m_minimizationStrategy;
		hair->m_params.m_rodParams.m_tolerance = it->second.m_minimizationTolerance;
		hair->m_params.m_rodParams.m_maxIter = it->second.m_minimizationMaxIter;

//		generate hair strands
		HairGenerator generator;
		generator.generateHair( *scene.m_renderObjects[idx], it->second.m_faceList, *hair );

		it->second.m_id = scene.m_hairs.size();
		scene.m_hairs.push_back( hair );
	}

	ifs.close();
	return true;
}


bool SceneLoader::loadMesh( const char* fileName, Mesh& mesh )
{
	ObjLoader loader;
	if ( !loader.loadFile( fileName ) )
	{
		std::cerr << "Error loading mesh from file " << fileName << std::endl;
		return false;
	}

	mesh.m_vertices.resize( loader.getVertices().size() );
	std::copy( loader.getVertices().begin(), loader.getVertices().end(), mesh.m_vertices.begin() );
	mesh.m_vindices.resize( loader.getVIndices().size() );
	std::copy( loader.getVIndices().begin(), loader.getVIndices().end(), mesh.m_vindices.begin() );

	if (    loader.getNormals().size() != loader.getVertices().size()
	     || loader.getVNIndices().size() != loader.getVIndices().size() )
	{
		Mesh::computeNormals( Mesh::ShadingMode::Gouraud, mesh );
	} else
	{
		mesh.m_normals.resize( loader.getNormals().size() );
		std::copy( loader.getNormals().begin(), loader.getNormals().end(), mesh.m_normals.begin() );
	}

	return true;
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
    // we require the standart "C" locale for numeric strings
    // otherwise parsing real numbers may fail depending on the machine setup
    char* savedLocale = std::setlocale(LC_NUMERIC, "C");

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

		success = false;
		std::cerr << "Unknown token[" << token << "] found at line "<< m_lineNumber << std::endl;
		break;
	}

    // restore locale
    setlocale(LC_NUMERIC, savedLocale);

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

		if (token.compare("length") == 0)
		{
			line.erase(0, token.length());
			object.m_length = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("lengthVariance") == 0)
		{
			line.erase(0, token.length());
			object.m_lengthVariance = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("helicalRadius") == 0)
		{
			line.erase(0, token.length());
			object.m_helicalRadius = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("helicalPitch") == 0)
		{
			line.erase(0, token.length());
			object.m_helicalPitch = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("density") == 0)
		{
			line.erase(0, token.length());
			object.m_density = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("thickness") == 0)
		{
			line.erase(0, token.length());
			object.m_thickness = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("nParticles") == 0)
		{
			line.erase(0, token.length());
			object.m_nParticles = BasicParser::parseInt(line);
			continue;
		}

		if (token.compare("netForce") == 0)
		{
			line.erase(0, token.length());
			BasicParser::parseVec3D(line, object.m_netForce);
			continue;
		}

		if (token.compare("drag") == 0)
		{
			line.erase(0, token.length());
			object.m_drag = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("resolveCollisions") == 0)
		{
			line.erase(0, token.length());
			object.m_resolveCollisions = BasicParser::parseInt(line);
			continue;
		}

		if (token.compare("resolveSelfInteractions") == 0)
		{
			line.erase(0, token.length());
			object.m_resolveSelfInterations = BasicParser::parseInt(line);
			continue;
		}

		if (token.compare("selfInteractionDist") == 0)
		{
			line.erase(0, token.length());
			object.m_selfInterationDist = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("selfStiction") == 0)
		{
			line.erase(0, token.length());
			object.m_selfStiction = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("selfRepulsion") == 0)
		{
			line.erase(0, token.length());
			object.m_selfRepulsion = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("pbdIter") == 0)
		{
			line.erase(0, token.length());
			object.m_pbdIter = BasicParser::parseInt(line);
			continue;
		}

		if (token.compare("bendStiffness") == 0)
		{
			line.erase(0, token.length());
			object.m_bendStiffness = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("twistStiffness") == 0)
		{
			line.erase(0, token.length());
			object.m_twistStiffness = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("maxElasticForce") == 0)
		{
			line.erase(0, token.length());
			object.m_maxElasticForce = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("minimizationMethod") == 0)
		{
			line.erase(0, token.length());
			object.m_minimizationStrategy = ElasticRodParams::MinimizationStrategy::None;

			std::string word = BasicParser::parseWord(BasicParser::ltrim(line));
			if (word.compare("newton") == 0)
			{
				object.m_minimizationStrategy =  ElasticRodParams::MinimizationStrategy::Newton;
			}
			else if (word.compare("bfgs") == 0)
			{
				object.m_minimizationStrategy =  ElasticRodParams::MinimizationStrategy::BFGS;
			}

			continue;
		}

		if (token.compare("minimizationTolerance") == 0)
		{
			line.erase(0, token.length());
			object.m_minimizationTolerance = BasicParser::parseDouble(line);
			continue;
		}

		if (token.compare("minimizationMaxIter") == 0)
		{
			line.erase(0, token.length());
			object.m_minimizationMaxIter = BasicParser::parseInt(line);
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
