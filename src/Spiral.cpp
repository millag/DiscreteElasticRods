#include "Spiral.h"

#include <QElapsedTimer>

Spiral::Spiral(): m_object(NULL)
{
    m_radius = 0.2;
    m_lenght = 2.0;

    m_bendStiffness = 10.0;
    m_twistStiffness = 10.0;
    m_maxForce = 5000;

    m_nParticles = 21;
    m_nIterations = 10;

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
    std::vector<bool> isFixed(m_nParticles, 0);

    mg::Vec3D center(m_object->getPosition().m_x, m_object->getPosition().m_y, m_object->getPosition().m_z);
    mg::Vec3D dirx = mg::Vec3D(1, 0, 0);
    mg::Vec3D diry = mg::Vec3D(0, -1, 0);
    mg::Vec3D dirz = mg::Vec3D(0, 0, 1);

    mg::Real step = 0.1;
    restpos[0] = center;
    restpos[1] = restpos[0] + step * diry;

    mg::Real angle = m_lenght / ((m_nParticles - 1) * m_radius);
    mg::Real len = step;
    for (unsigned i = 2; i < restpos.size(); ++i)
    {
        restpos[i] = std::cos((i - 1) * angle) * dirx * m_radius + std::sin((i - 1) * angle) * dirz * m_radius + center - dirx * m_radius;
        restpos[i] += (i - 1) * step * diry;
        len += mg::length(restpos[i] - restpos[i - 1]);
    }
    isFixed[0] = 1;
//    isFixed[m_nParticles - 1] = 1;


    mg::Vec3D start = restpos[0];
    mg::Vec3D end = start + len * diry;
    mg::Real t = 0.0;
    for (unsigned i = 0; i < pos.size(); ++i)
    {
        t = (mg::Real)(i) / (m_nParticles - 1);
        pos[i] = (1 - t) * start + t * end;

        vel[i].zero();
        mass[i] = 10;
    }

    for (unsigned i = 0; i < 1; ++i)
    {
        ElasticRod* strand = new ElasticRod(m_bendStiffness, m_twistStiffness, m_nIterations, m_maxForce, ElasticRod::BFGS_NUMERIC);
        strand->init(restpos, dirx, restpos, vel, mass, twistAngle, isFixed);
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
        (*it)->m_ppos[0].set(m_object->getPosition().m_x, m_object->getPosition().m_y, m_object->getPosition().m_z);
        (*it)->update(dt);
    }

    std::cout << "TIME update ms: " << chronometer.elapsed() << std::endl;
    chronometer.restart();
}
