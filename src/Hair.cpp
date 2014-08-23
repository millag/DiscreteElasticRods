#include "Hair.h"
#include <QElapsedTimer>
#include "Utils.h"

HairParams::HairParams(): m_rodParams(NULL)
{ }

HairParams::~HairParams()
{
    if (m_rodParams != NULL)
    {
        delete m_rodParams;
    }
}

Hair::Hair():m_object(NULL)
{
//    init default hair params
    m_params = new HairParams();
    m_params->m_length = 4;
    m_params->m_lengthVariance = 0;
    m_params->m_helicalRadius = 0.2;
    m_params->m_helicalPitch = 0.1;
    m_params->m_density = 0.5;
    m_params->m_thickness = 0.07;
    m_params->m_nParticles = 25;

    m_params->m_gravity.set(0, -9.81, 0);
    m_params->m_drag = 0.001;

    mg::Real bendStiffness = 0.01;
    mg::Real twistStiffness = 0.01;
    unsigned pbdIterations = 5;
    mg::Real maxElasticForce = 1000;
    m_params->m_rodParams = new RodParams(bendStiffness, twistStiffness, pbdIterations, maxElasticForce);
}


Hair::~Hair()
{
    reset();
    delete m_params;
}

void Hair::reset()
{
//    delete hair strands
    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        delete (*it);
    }
    m_strands.clear();
    m_findices.clear();
    m_vindices.clear();
    m_object = NULL;
}

void Hair::update(mg::Real dt)
{
    assert(m_object != NULL);

    const Mesh* mesh = m_object->getMesh();

    QElapsedTimer chronometer;
    chronometer.start();

#ifdef MULTI_THREAD
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < m_vindices.size(); ++i)
    {
        m_strands[i]->m_ppos[0] = mg::transform_point(m_object->getTransform(), mesh->m_vertices[ m_vindices[i] ]);
        updateRod(*m_strands[i], dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}

/* semi-implicit Euler with Verlet scheme for velocity update */
void Hair::updateRod(ElasticRod& rod, mg::Real dt) const
{
    std::vector<mg::Vec3D> forces(rod.m_ppos.size(), mg::Vec3D(0,0,0));
    rod.accumulateInternalElasticForces(forces);
    accumulateExternalForces(rod, forces);

//    integrate centerline - semi- implicit Euler
    std::vector<mg::Vec3D> prevPos(rod.m_ppos.size());
    for (unsigned i = 0; i < rod.m_ppos.size(); ++i)
    {
        prevPos[i] = rod.m_ppos[i];
        rod.m_pvel[i] += dt * forces[i] / rod.m_pmass[i];
        rod.m_ppos[i] += rod.m_pvel[i] * dt;
    }

    rod.enforceInternalConstraints();

//    update velocities:
    for (unsigned i = 0; i < rod.m_ppos.size(); ++i)
    {
        rod.m_pvel[i] = (rod.m_ppos[i] - prevPos[i]) / dt;
    }

    rod.updateCurrentState();
}

void Hair::accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
    for (unsigned i = 0; i < o_forces.size(); ++i)
    {
        if (!rod.m_isClamped.count(i))
        {
            o_forces[i] += m_params->m_gravity * rod.m_pmass[i] - m_params->m_drag * rod.m_pvel[i];
        }
    }
}
