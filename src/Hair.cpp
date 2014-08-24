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

Hair::Hair():m_object(NULL), m_grid(NULL)
{
//    init default hair params
    m_params = new HairParams();
    m_params->m_length = 3;
    m_params->m_lengthVariance = 0.9;
    m_params->m_helicalRadius = 0.2;
    m_params->m_helicalPitch = 0.1;
    m_params->m_density = 0.09;
    m_params->m_thickness = 0.07;
    m_params->m_nParticles = 12;

    m_params->m_coulombFriction = 0.2;
    m_params->m_selfCollisionDist = 0.6;
    m_params->m_selfFriction = 0.03;
    m_params->m_selfRepulsion = 0.004;

    m_params->m_gravity.set(0, -9.81, 0);
    m_params->m_drag = 0.0005;

    m_params->m_resolveSelfCollisions = 1;
    m_params->m_pbdIter = 4;

    mg::Real bendStiffness = 0.003;
    mg::Real twistStiffness = 0.0001;
    mg::Real maxElasticForce = 1000;
    m_params->m_rodParams = new RodParams(bendStiffness, twistStiffness, maxElasticForce);
}


Hair::~Hair()
{
    reset();
    delete m_params;
}

void Hair::reset()
{
    if (m_grid != NULL)
    {
        delete m_grid;
    }
    typedef std::vector<ElasticRod*>::iterator Iter;
    for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        delete (*it);
    }
    m_strands.clear();
    m_findices.clear();
    m_vindices.clear();
    m_object = NULL;
    m_grid = NULL;
}

void Hair::createGrid()
{
    mg::Vec3D offset = (m_params->m_length + m_params->m_lengthVariance) * mg::Vec3D(1, 1, 1);
    m_volume.reshape(m_object->getAABB().getVMin() - offset, m_object->getAABB().getVMax() + offset);

    m_grid = new VoxelGrid(m_volume, m_volume.getWidth() / m_params->m_selfCollisionDist);
    m_grid->initialize();
}

void Hair::resetGrid()
{
    m_grid->reset();

    typedef std::vector<ElasticRod*>::iterator Iter;
    for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        ElasticRod* rod = *it;
        for (unsigned i = 1; i < rod->m_ppos.size(); ++i)
        {
            m_grid->insertDensity(rod->m_ppos[i], 1);
        }
    }

    for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        ElasticRod* rod = *it;
        for (unsigned i = 1; i < rod->m_ppos.size(); ++i)
        {
            m_grid->insertVelocity(rod->m_ppos[i], rod->m_pvel[i]);
        }
    }

}

void Hair::update(mg::Real dt)
{
    assert(m_object != NULL);

    QElapsedTimer chronometer;
    chronometer.start();

    if (m_params->m_resolveSelfCollisions)
    {
        if (m_grid == NULL)
        {
            createGrid();
        } else
        {
            resetGrid();
        }
    }

    const Mesh* mesh = m_object->getMesh();


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
    for (unsigned i = 1; i < rod.m_ppos.size(); ++i)
    {
        prevPos[i] = rod.m_ppos[i];
        rod.m_pvel[i] += dt * forces[i] / rod.m_pmass[i];
        rod.m_ppos[i] += rod.m_pvel[i] * dt;
    }

    for (unsigned k = 0; k < m_params->m_pbdIter; ++k)
    {
        applyCollisionConstraintsIteration(rod);
        rod.applyInternalConstraintsIteration();
    }

//    velocity correction using Verlet scheme:
    for (unsigned i = 1; i < rod.m_ppos.size(); ++i)
    {
        rod.m_pvel[i] = (rod.m_ppos[i] - prevPos[i]) / dt;
    }

    rod.updateCurrentState();
}

//void Hair::accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
//{
//    for (unsigned i = 1; i < o_forces.size(); ++i)
//    {
//        if (!rod.m_isClamped.count(i))
//        {
//            o_forces[i] += m_params->m_gravity * rod.m_pmass[i] - m_params->m_drag * rod.m_pvel[i];
//        }
//    }
//}

void Hair::applyCollisionConstraintsIteration(ElasticRod& rod) const
{
    mg::Vec3D collision_p, normal;
    for ( unsigned i = 1; i < rod.m_ppos.size(); ++i)
    {
        mg::Real dist = m_object->distanceToSurface(rod.m_ppos[i], collision_p, normal);
        if (dist < m_params->m_thickness)
        {
            rod.m_ppos[i] = collision_p + m_params->m_thickness * normal;
        }

    }
}

void Hair::accumulateExternalForces(ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
    mg::Vec3D p1(0,0,0), p2(0,0,0), pressureGradient(0,0,0);
    mg::Real dr = m_params->m_selfCollisionDist;
    mg::Real density_p1, density_p2;
    for (unsigned i = 1; i < o_forces.size(); ++i)
    {
//        gravity
        o_forces[i] += m_params->m_gravity * rod.m_pmass[i];
//        drag
        o_forces[i] += -m_params->m_drag * rod.m_pvel[i];


//        self interactions
        if (!m_params->m_resolveSelfCollisions)
        {
            continue;
        }

//            calculate self friction
        m_grid->getInterpolatedDensity(rod.m_ppos[i], density_p1);
        if ( density_p1 > mg::ERR)
        {
            m_grid->getInterpolatedVelocity(rod.m_ppos[i], p1);
            rod.m_pvel[i] =  (1 - m_params->m_selfFriction) * rod.m_pvel[i] + m_params->m_selfFriction * p1 / density_p1;
        }

//        self repulsion
        p1.set(rod.m_ppos[i][0] - dr, rod.m_ppos[i][1], rod.m_ppos[i][2]);
        p2.set(rod.m_ppos[i][0] + dr, rod.m_ppos[i][1], rod.m_ppos[i][2]);
        m_grid->getInterpolatedDensity(p1, density_p1);
        m_grid->getInterpolatedDensity(p2, density_p2);
        pressureGradient[0] = ( density_p1 - density_p2 ) * 0.5 / dr;

        p1.set(rod.m_ppos[i][0], rod.m_ppos[i][1] - dr, rod.m_ppos[i][2]);
        p2.set(rod.m_ppos[i][0], rod.m_ppos[i][1] + dr, rod.m_ppos[i][2]);
        m_grid->getInterpolatedDensity(p1, density_p1);
        m_grid->getInterpolatedDensity(p2, density_p2);
        pressureGradient[1] = ( density_p1 - density_p2 ) * 0.5 / dr;

        p1.set(rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] - dr);
        p2.set(rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] + dr);
        m_grid->getInterpolatedDensity(p1, density_p1);
        m_grid->getInterpolatedDensity(p2, density_p2);
        pressureGradient[2] = ( density_p1 - density_p2) * 0.5 / dr;

        o_forces[i] += m_params->m_selfRepulsion * pressureGradient;

    }
}

//for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
//{
//    Strand& strand = *it;
//    for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
//    {
////            calculate self friction
//        strand.m_pvel[i] =  (1 - m_selfFriction) * strand.m_pvel[i] +
//                            m_selfFriction * m_grid->getInterpolatedVelocity(strand.m_ppos[i]) / m_grid->getInterpolatedDensity(strand.m_ppos[i]);

////            calculate self repulsion
//        p1.set(strand.m_ppos[i][0] - dr, strand.m_ppos[i][1], strand.m_ppos[i][2]);
//        p2.set(strand.m_ppos[i][0] + dr, strand.m_ppos[i][1], strand.m_ppos[i][2]);
//        pressureGradient[0] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        p1.set(strand.m_ppos[i][0], strand.m_ppos[i][1] - dr, strand.m_ppos[i][2]);
//        p2.set(strand.m_ppos[i][0], strand.m_ppos[i][1] + dr, strand.m_ppos[i][2]);
//        pressureGradient[1] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        p1.set(strand.m_ppos[i][0], strand.m_ppos[i][1], strand.m_ppos[i][2] - dr);
//        p2.set(strand.m_ppos[i][0], strand.m_ppos[i][1], strand.m_ppos[i][2] + dr);
//        pressureGradient[2] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        strand.m_pvel[i] = 0.5 * (2 * strand.m_pvel[i] +  m_selfRepulsion * dt * pressureGradient / strand.m_pmass[i]);
//    }
//}

//        add friction
//        if (rod.m_isClamped.count(i))
//        {
//            continue;
//        }
//        mg::Vec3D collision_p, normal, vt, vn;
//        mg::Real dist = m_object->distanceToSurface(rod.m_ppos[i], collision_p, normal);
//        if (dist > m_params->m_thickness + 0.0001)
//        {
//            continue;
//        }

//        vn = mg::dot(normal, rod.m_pvel[i]) * normal;
//        vt = rod.m_pvel[i] + vn;

//        if (vt.length_squared() > mg::ERR)
//        {
//            rod.m_pvel[i] -= std::max(1 - m_params->m_coulombFriction * vn.length() / vt.length(), (mg::Real)0) * vt;
//        }



//// =========================== long and wavy ========================================
//    m_params = new HairParams();
//    m_params->m_length = 6;
//    m_params->m_lengthVariance = 1;
//    m_params->m_helicalRadius = 0.2;
//    m_params->m_helicalPitch = 0.13;
//    m_params->m_density = 0.01;
//    m_params->m_thickness = 0.07;
//    m_params->m_nParticles = 25;

//    m_params->m_coulombFriction = 0.2;
//    m_params->m_selfCollisionDist = 0.5;
//    m_params->m_selfFriction = 0.001;
//    m_params->m_selfRepulsion = 0.0001;

//    m_params->m_gravity.set(0, -9.81, 0);
//    m_params->m_drag = 0.0001;

//    m_params->m_resolveSelfCollisions = 0;
//    m_params->m_pbdIter = 4;

//    mg::Real bendStiffness = 0.0003;
//    mg::Real twistStiffness = 0.0001;
//    mg::Real maxElasticForce = 1000;
//    m_params->m_rodParams = new RodParams(bendStiffness, twistStiffness, maxElasticForce);
