#include "HairStrand.h"

#include <QElapsedTimer>

HairStrand::HairStrand():m_object(NULL)
{
    m_length = 4.0;
    m_lengthVariance = 1.0;
    m_density = 0.001;
    m_thickness = 0.001;
    m_bendStiffness = 5.0;
    m_twistStiffness = 2.0;
    m_maxForce = 500;

    m_nParticles = 15;
    m_nIterations = 4;

    m_strands.reserve(100);
}

HairStrand::~HairStrand()
{
//    delete hair strands
    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        delete (*it);
    }
}

void HairStrand::initialize(const RenderObject* object)
{
    assert(object != NULL);
    m_object = object;

    mg::Vec3D start(m_object->getPosition().m_x, m_object->getPosition().m_y, m_object->getPosition().m_z);
    mg::Vec3D end = start + mg::Vec3D(1, 0, 0) * m_length;

    std::vector<mg::Vec3D> restpos(m_nParticles);
    std::vector<mg::Vec3D> pos(m_nParticles);
    std::vector<mg::Vec3D> vel(m_nParticles);
    std::vector<mg::Real> mass(m_nParticles);
    std::vector<double> twistAngle(m_nParticles - 1, 0.0);
    std::set<unsigned> isClamped;

    mg::Real t = 0.0;
    for (unsigned i = 0; i < pos.size(); ++i)
    {
        t = (mg::Real)(i) / (m_nParticles - 1);
        pos[i] = (1 - t) * start + t * end;

        vel[i].zero();
        mass[i] = 1;
    }

    mg::Real step = mg::length( end - start ) / (m_nParticles - 1);
    restpos[0] = start;
    restpos[1] = restpos[0] + step * mg::Vec3D(1, 0, 0);
    mg::Vec3D dir = mg::normalize(mg::Vec3D(1, 1, 1));
    for (unsigned i = 2; i < restpos.size(); ++i)
    {
        if ((i - 1) % 2)
        {
            dir[1] = -dir[1];
        }
        if (i % 2)
        {
            dir[2] = -dir[2];
        }
        restpos[i] = restpos[i - 1] + step * dir;
    }
    isClamped.insert(0);
//    isClamped.insert(m_nParticles - 1);


    for (unsigned i = 0; i < 1; ++i)
    {
        ElasticRod* strand = new ElasticRod(m_bendStiffness, m_twistStiffness, m_nIterations, m_maxForce, ElasticRod::NEWTON);
        strand->init(restpos, mg::Vec3D(0,1,0), restpos, vel, mass, twistAngle, isClamped);
        m_strands.push_back(strand);
    }
}

void HairStrand::update(mg::Real dt)
{
    QElapsedTimer chronometer;
    chronometer.start();

    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        (*it)->m_ppos[0].set(m_object->getPosition().m_x, m_object->getPosition().m_y, m_object->getPosition().m_z);
        (*it)->update(dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}
