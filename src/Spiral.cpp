#include "Spiral.h"

#include <QElapsedTimer>

Spiral::Spiral(): m_object(NULL)
{
    m_radius = 0.2;
    m_lenght = 4.0;

    m_bendStiffness = 0.1;
    m_twistStiffness = 0.01;
    m_maxForce = 10;

    m_nParticles = 25;
    m_nIterations = 4;

    m_strands.reserve(100);
}

Spiral::~Spiral()
{
    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        delete (*it);
    }
}

void Spiral::init(const RenderObject* object)
{
    assert(object != NULL);
    m_object = object;

    std::vector<mg::Vec3D> restpos(m_nParticles);
    std::vector<mg::Vec3D> pos(m_nParticles);
    std::vector<mg::Vec3D> vel(m_nParticles);
    std::vector<mg::Real> mass(m_nParticles);
    std::vector<double> twistAngle(m_nParticles - 1, 0.0);
    std::set<unsigned> isClamped;

    mg::Vec3D center(m_object->getPosition());
    mg::Vec3D dirx = mg::Vec3D(1, 0, 0);
    mg::Vec3D diry = mg::Vec3D(0, -1, 0);
    mg::Vec3D dirz = mg::Vec3D(0, 0, 1);


    mg::Real angle = m_lenght / ((m_nParticles - 1) * m_radius);
    mg::Real step = 0.1;
    mg::Real len = 0;
    for (unsigned i = 0; i < restpos.size(); ++i)
    {
        restpos[i] = (std::cos(i * angle) * dirx  + std::sin(i * angle) * dirz - dirx) * m_radius + center;
        restpos[i] += diry * i * step;
        if (i > 0)
        {
            len += mg::length(restpos[i] - restpos[i - 1]);
        }
    }
    isClamped.insert(0);
//    isClamped.insert(m_nParticles - 1);

    mg::Vec3D u0 = mg::cross((restpos[1] - restpos[0]), dirx);
    u0 = mg::normalize( mg::cross(u0, (restpos[1] - restpos[0])) );

    mg::Vec3D start = restpos[0];
    mg::Vec3D end = start + len * dirx;
    mg::Real t = 0.0;
    for (unsigned i = 0; i < pos.size(); ++i)
    {
        t = (mg::Real)(i) / (m_nParticles - 1);
        pos[i] = (1 - t) * start + t * end;

        vel[i].zero();
        mass[i] = 0.01;
    }

    for (unsigned i = 0; i < 1; ++i)
    {
        ElasticRod* strand = new ElasticRod(m_bendStiffness, m_twistStiffness, m_nIterations, m_maxForce, ElasticRod::BFGS);
        strand->init(restpos, u0, restpos, vel, mass, twistAngle, isClamped);
        m_strands.push_back(strand);
    }
}

void Spiral::update(mg::Real dt)
{
    QElapsedTimer chronometer;
    chronometer.start();

    typedef std::vector<ElasticRod*>::const_iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        (*it)->m_ppos[0] = m_object->getPosition();
        (*it)->update(dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}
