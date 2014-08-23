#include "Hair.h"
#include <QElapsedTimer>
#include "config.h"

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
        unsigned idx = m_vindices[i];
        m_strands[i]->m_ppos[0] = mg::transform_point(m_object->getTransform(), mesh->m_vertices[idx]);
        m_strands[i]->update(dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}
