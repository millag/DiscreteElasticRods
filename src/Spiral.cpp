#include "Spiral.h"
#include "Utils.h"

Spiral::Spiral(): m_object(NULL)
{
    m_radius = 0.2;
    m_lenght = 4.0;

    m_bendStiffness = 0.4;
    m_twistStiffness = 0.7;
    m_maxForce = 1000;

    m_offset = 3;

    m_nParticles = 20;
    m_pbdIter = 4;

    m_rodParams1 = new ElasticRodParams(m_bendStiffness, m_twistStiffness, m_maxForce, ElasticRodParams::NEWTON);
    m_rodParams2 = new ElasticRodParams(m_bendStiffness, m_twistStiffness, m_maxForce, ElasticRodParams::BFGS);
    m_rodParams3 = new ElasticRodParams(m_bendStiffness, m_twistStiffness, m_maxForce, ElasticRodParams::NONE);
    m_strands.reserve(10);
}

Spiral::~Spiral()
{
    typedef std::vector<ElasticRod*>::iterator Iter;
    for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        delete (*it);
    }
    m_strands.clear();
    delete m_rodParams1;
    delete m_rodParams2;
    delete m_rodParams3;
}

void Spiral::init(const RenderObject* object)
{
    assert(object != NULL);
    m_object = object;

    std::vector<mg::Vec3D> restpos(m_nParticles);
    std::vector<mg::Vec3D> pos(m_nParticles);
    std::vector<mg::Vec3D> vel(m_nParticles);
    std::vector<mg::Real> mass(m_nParticles);
    ColumnVector theta(m_nParticles - 1);
    std::set<unsigned> isClamped;

    mg::Vec3D center(m_object->getCenter());
    mg::Vec3D dirx(1, 0, 0);
    mg::Vec3D diry(0, -1, 0);
    mg::Vec3D dirz(0, 0, 1);


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
        mass[i] = 0.05;
        if (i < (unsigned)theta.size())
        {
            theta(i) = 0;
        }
    }

    ElasticRod* rod = new ElasticRod(m_rodParams1);
//        strand->init(pos, diry, pos, vel, mass, theta, isClamped);
    rod->init(restpos, u0, pos, vel, mass, theta, isClamped);
    m_strands.push_back(rod);

    for (unsigned i = 0; i < pos.size(); ++i)
    {
        pos[i] -= m_offset * dirx;
    }

    rod = new ElasticRod(m_rodParams2);
//        strand->init(pos, diry, pos, vel, mass, theta, isClamped);
    rod->init(restpos, u0, pos, vel, mass, theta, isClamped);
    m_strands.push_back(rod);


    for (unsigned i = 0; i < pos.size(); ++i)
    {
        pos[i] += 2 * m_offset * dirx;
    }

    rod = new ElasticRod(m_rodParams3);
//        strand->init(pos, diry, pos, vel, mass, theta, isClamped);
    rod->init(restpos, u0, pos, vel, mass, theta, isClamped);
    m_strands.push_back(rod);
}

void Spiral::update(mg::Real dt)
{
    mg::Vec3D center = m_object->getCenter();

    m_strands[0]->m_ppos[0] = center;
    updateRod(*m_strands[0], dt);

    m_strands[1]->m_ppos[0] = center - m_offset * mg::Vec3D(1, 0, 0);
    updateRod(*m_strands[1], dt);

    m_strands[2]->m_ppos[0] = center + m_offset * mg::Vec3D(1, 0, 0);
    updateRod(*m_strands[2], dt);
}

/* semi-implicit Euler with Verlet scheme for velocity update */
void Spiral::updateRod(ElasticRod& rod, mg::Real dt) const
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

    for (unsigned k = 0; k < m_pbdIter; ++k)
    {
        rod.applyInternalConstraintsIteration();
    }

//    velocity correction using Verlet scheme:
    for (unsigned i = 0; i < rod.m_ppos.size(); ++i)
    {
        rod.m_pvel[i] = (rod.m_ppos[i] - prevPos[i]) / dt;
    }

    rod.updateCurrentState();
}

void Spiral::accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
    for (unsigned i = 0; i < o_forces.size(); ++i)
    {
        if (!rod.m_isClamped.count(i))
        {
            o_forces[i] += mg::GRAVITY * rod.m_pmass[i] - 0.001 * rod.m_pvel[i];
        }
    }
}


//void Spiral::init(const RenderObject* object)
//{
//    assert(object != NULL);
//    m_object = object;

//    std::vector<mg::Vec3D> restpos(m_nParticles);
//    std::vector<mg::Vec3D> pos(m_nParticles);
//    std::vector<mg::Vec3D> vel(m_nParticles);
//    std::vector<mg::Real> mass(m_nParticles);
//    ColumnVector theta(m_nParticles - 1);
//    std::set<unsigned> isClamped;

//    mg::Vec3D center(m_object->getCenter());
//    mg::Vec3D dirx = mg::Vec3D(1, 0, 0);
//    mg::Vec3D diry = mg::Vec3D(0, 1, 0);
//    mg::Vec3D dirz = mg::Vec3D(0, 0, 1);

//    mg::Real radius = m_lenght / mg::Constants::pi();
//    mg::Real angle = mg::Constants::pi() / (m_nParticles - 1);
//    for (unsigned i = 0; i < restpos.size(); ++i)
//    {
//        unsigned j = i;//(m_nParticles - 1) - i;
//        pos[i] = (std::cos(j * angle) * dirx  + std::sin(j * angle) * diry - dirx) * radius + center;
//    }
//    isClamped.insert(0);
////    isClamped.insert(m_nParticles - 1);

//    mg::Vec3D start = pos[0];
//    mg::Vec3D end = start + m_lenght * dirx;
//    mg::Real t = 0.0;
//    for (unsigned i = 0; i < pos.size(); ++i)
//    {
//        t = (mg::Real)(i) / (m_nParticles - 1);
//        restpos[i] = (1 - t) * start + t * end;

//        vel[i].zero();
//        mass[i] = 0.001;
//        if (i < (unsigned)theta.size())
//        {
//            theta(i) = 0;
//        }
//    }

//    for (unsigned i = 0; i < 1; ++i)
//    {
//        ElasticRod* strand = new ElasticRod(m_rodParams);
////        strand->init(pos, diry, pos, vel, mass, theta, isClamped);
//        strand->init(restpos, -dirz, pos, vel, mass, theta, isClamped);
//        m_strands.push_back(strand);
//    }
//}
