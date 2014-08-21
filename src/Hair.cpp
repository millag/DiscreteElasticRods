#include "Hair.h"
#include <QElapsedTimer>

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
    m_params = new HairParams();
    m_params->m_length = 4;
    m_params->m_lengthVariance = 1;
    m_params->m_helicalRadius = 0.2;
    m_params->m_helicalPitch = 0.1;
    m_params->m_density = 0.001;
    m_params->m_thickness = 0.1;
    m_params->m_nParticles = 15;

    mg::Real bendStiffness = 5.0;
    mg::Real twistStiffness = 2.0;
    unsigned pbdIterations = 4;
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
    m_object = NULL;
}

void Hair::update(mg::Real dt)
{
    assert(m_object != NULL);

    QElapsedTimer chronometer;
    chronometer.start();

    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
//        TODO: need to update pos via transform
        (*it)->m_ppos[0] = m_object->getPosition();
        (*it)->update(dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}
