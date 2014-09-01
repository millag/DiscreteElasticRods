#include "HairGenerator.h"
#include <QElapsedTimer>
#include "Utils.h"
#include "config.h"


HairGenerator::HairGenerator()
{ }

HairGenerator::~HairGenerator()
{ }

void HairGenerator::generateCurlyHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair)
{
    assert(object != NULL);
    assert(object->getMesh() != NULL);

    const Mesh* mesh = object->getMesh();

    o_hair.reset();
    o_hair.m_object = object;
    o_hair.m_findices = findices;
    o_hair.m_vindices.reserve( findices.size() );

    std::pair<std::set<unsigned>::const_iterator, bool> hasIdx;
    std::set<unsigned> vIdxSet;
    for (unsigned i = 0; i < findices.size(); ++i)
    {
        unsigned idx = mesh->getPrimitiveOffset(findices[i]);
        for (unsigned j = 0; j < mesh->getNVerticesPerPrimitive(); ++j)
        {
            hasIdx = vIdxSet.insert(mesh->m_vindices[idx + j]);
            if (hasIdx.second)
            {
                o_hair.m_vindices.push_back( *hasIdx.first );
            }
        }
    }

    o_hair.m_strands.resize( o_hair.m_vindices.size() );

#ifdef MULTI_THREAD
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < o_hair.m_vindices.size(); ++i)
    {
        unsigned idx = o_hair.m_vindices[i];
        mg::Vec3D p = mg::transform_point(object->getTransform(), mesh->m_vertices[ idx ]);
        mg::Vec3D n = mg::transform_vector(object->getTransform(), mesh->m_normals[ idx ]);
        mg::Vec3D u = mg::EY;
        if (fabs(1 - fabs(mg::dot(n, u))) < mg::ERR)
        {
            u = mg::EX;
        }

        ElasticRod* rod = new ElasticRod(o_hair.m_params->m_rodParams);
        generateHelicalRod(*o_hair.m_params, p, n, u, *rod);
        o_hair.m_strands[i] = rod;
    }

    o_hair.initialize();
}

void HairGenerator::generateStraightHair(const RenderObject* object, const std::vector<unsigned>& findices, Hair& o_hair)
{
    assert(object != NULL);
    assert(object->getMesh() != NULL);

    const Mesh* mesh = object->getMesh();

    o_hair.reset();
    o_hair.m_object = object;
    o_hair.m_findices = findices;
    o_hair.m_vindices.reserve( findices.size() );

    std::pair<std::set<unsigned>::const_iterator, bool> hasIdx;
    std::set<unsigned> vIdxSet;
    for (unsigned i = 0; i < findices.size(); ++i)
    {
        unsigned idx = mesh->getPrimitiveOffset(findices[i]);
        for (unsigned j = 0; j < mesh->getNVerticesPerPrimitive(); ++j)
        {
            hasIdx = vIdxSet.insert(mesh->m_vindices[idx + j]);
            if (hasIdx.second)
            {
                o_hair.m_vindices.push_back( *hasIdx.first );
            }
        }
    }

    o_hair.m_strands.resize( o_hair.m_vindices.size() );

#ifdef MULTI_THREAD
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < o_hair.m_vindices.size(); ++i)
    {
        unsigned idx = o_hair.m_vindices[i];
        mg::Vec3D p = mg::transform_point(object->getTransform(), mesh->m_vertices[ idx ]);
        mg::Vec3D n = mg::transform_vector(object->getTransform(), mesh->m_normals[ idx ]);
        mg::Vec3D u = mg::EY;
        if (fabs(1 - fabs(mg::dot(n, u))) < mg::ERR)
        {
            u = mg::EX;
        }

        ElasticRod* rod = new ElasticRod(o_hair.m_params->m_rodParams);
        generateStraightRod(*o_hair.m_params, p, n, u, *rod);
        o_hair.m_strands[i] = rod;
    }

    o_hair.initialize();
}

void HairGenerator::generateHelicalRod(const HairParams& params,
                                    const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                                    ElasticRod& o_rod)
{
    o_rod.m_ppos.resize(params.m_nParticles);
    o_rod.m_pvel.resize(params.m_nParticles);
    o_rod.m_pmass.resize(params.m_nParticles);
    o_rod.m_theta.set_size(params.m_nParticles - 1);
    o_rod.m_isClamped.clear();

    mg::Vec3D dirv = mg::normalize(mg::cross(n, u));
    mg::Vec3D dirn = mg::normalize(n);
    mg::Vec3D diru = mg::normalize(mg::cross(dirv, dirn));

    mg::Real length = params.m_length + mg::randf(-params.m_lengthVariance, params.m_lengthVariance);
    mg::Real volume = length * params.m_thickness * params.m_thickness * mg::Constants::pi();
    mg::Real pmass = params.m_density * volume / (params.m_nParticles - 1);

    mg::Real sign = mg::randf(0.0, 1.0);
    sign =  (sign < 0.5)? -1.0 : 1.0;
    mg::Real angle = params.m_length / std::sqrt(params.m_helicalRadius * params.m_helicalRadius + params.m_helicalPitch * params.m_helicalPitch);
    angle /= (params.m_nParticles - 1);
    mg::Real phase = mg::randf(0, mg::Constants::two_pi());
    for (unsigned i = 0; i < o_rod.m_ppos.size(); ++i)
    {
//        calc point on unit circle
        o_rod.m_ppos[i] = (std::cos(phase + i *sign * angle) * diru  + std::sin(phase + i * sign * angle) * dirv - (std::cos(phase) * diru  + std::sin(phase) * dirv));
//        place point on helix
        o_rod.m_ppos[i] = o_rod.m_ppos[i]* params.m_helicalRadius + params.m_helicalPitch * (i * angle) * dirn;
//        move helix to position
        o_rod.m_ppos[i] += p;
//        init velocity and mass
        o_rod.m_pvel[i].zero();
        o_rod.m_pmass[i] =  pmass;
        if (i < static_cast<unsigned>(o_rod.m_theta.size()))
        {
            o_rod.m_theta(i) = 0;
        }
    }
    o_rod.m_isClamped.insert(0);
//    o_rod.m_isClamped.insert(params->m_nParticles - 1);

    mg::Vec3D e0 = o_rod.m_ppos[1] - o_rod.m_ppos[0];
    o_rod.m_u0 = mg::cross(e0, u);
    o_rod.m_u0 = mg::normalize( mg::cross(o_rod.m_u0, e0) );

    o_rod.init(o_rod.m_ppos, o_rod.m_u0, o_rod.m_ppos, o_rod.m_pvel, o_rod.m_pmass, o_rod.m_theta, o_rod.m_isClamped);
}


void HairGenerator::generateStraightRod(const HairParams& params,
                                     const mg::Vec3D& p, const mg::Vec3D& n, const mg::Vec3D& u,
                                     ElasticRod& o_rod)
{
    o_rod.m_ppos.resize(params.m_nParticles);
    o_rod.m_pvel.resize(params.m_nParticles);
    o_rod.m_pmass.resize(params.m_nParticles);
    o_rod.m_theta.set_size(params.m_nParticles - 1);
    o_rod.m_isClamped.clear();

    mg::Real length = params.m_length + mg::random_real(-params.m_lengthVariance, params.m_lengthVariance);
    mg::Real volume = length * params.m_thickness * params.m_thickness * mg::Constants::pi();
    mg::Real pmass = params.m_density * volume / (params.m_nParticles - 1);

//    generate straight line
    mg::Vec3D end = p + length * mg::normalize(n);
    mg::Real t = 0.0;
    for (unsigned i = 0; i < o_rod.m_ppos.size(); ++i)
    {
        t = (mg::Real)(i) / (params.m_nParticles - 1);
        o_rod.m_ppos[i] = (1 - t) * p + t * end;
//        init velocity and mass
        o_rod.m_pvel[i].zero();
        o_rod.m_pmass[i] = pmass;
        if (i < static_cast<unsigned>(o_rod.m_theta.size()))
        {
            o_rod.m_theta(i) = 0;
        }
    }
    o_rod.m_isClamped.insert(0);
//    o_rod.m_isClamped.insert(params->m_nParticles - 1);

    mg::Vec3D e0 = o_rod.m_ppos[1] - o_rod.m_ppos[0];
    o_rod.m_u0 = mg::cross(e0, u);
    o_rod.m_u0 = mg::normalize( mg::cross(o_rod.m_u0, e0) );

    o_rod.init(o_rod.m_ppos, o_rod.m_u0, o_rod.m_ppos, o_rod.m_pvel, o_rod.m_pmass, o_rod.m_theta, o_rod.m_isClamped);
}
