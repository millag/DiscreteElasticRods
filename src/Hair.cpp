#include "Hair.h"
#include <cmath>
#include <ngl/Util.h>
#include "Utils.h"

Hair::Hair(IServant* servant):
    m_servant(servant), m_grid(NULL)
{
    assert(m_servant != NULL);

    m_useDFTL = true;

//    long hair preset
    m_length = 2.5;
    m_lengthVariance = 0.8;
    m_density = 0.001;
    m_pmassFalloff = 0.0;
    m_thickness = 0.001;
    m_stiffness = 1;
    m_damping = 0.8;
    m_collisionFriction = 0.5;
//    air drag force is computed by the formula: dragCoefficient * 1/2 * airDensity * v^2 * areaCrossSecction
    m_airDrag = 10000 * 0.5 * 1.225 * 0.25 * ngl::PI * m_thickness * m_thickness;
    m_nParticles = 8;
    m_nIterations = 4;
    m_selfCollisionDist = 0.2;
    m_selfFriction = 0.4;
    m_selfRepulsion = 0.05;

//    short fur
//    m_length = 0.6;
//    m_lengthVariance = 0.4;
//    m_density = 0.001;
//    m_pmassFalloff = 0.1;
//    m_thickness = 0.001;
//    m_stiffness = 1;
//    m_damping = 1;
//    m_collisionFriction = 1;
////    air drag force is computed by the formula: dragCoefficient * 1/2 * airDensity * v^2 * areaCrossSecction
//    m_airDrag = 10000 * 0.5 * 1.225 * 0.25 * ngl::PI * m_thickness * m_thickness;
//    m_nParticles = 4;
//    m_nIterations = 4;
//    m_selfCollisionDist = 0.3;
//    m_selfFriction = 0.8;
//    m_selfRepulsion = 0.1;
}

Hair::~Hair()
{
    if (m_grid != NULL)
        delete m_grid;
}

//    TODO: implement Poison disc on selected faces
void Hair::initialize(const RenderObject& object, unsigned nMaxStrands)
{
    assert(nMaxStrands > 0);
    assert(m_nParticles > 1);

    m_strands.resize(nMaxStrands);

    const ngl::Mat4& t = object.getTransform();
    ngl::Vec4 initPos;
    ngl::Vec4 dir;
    typedef std::vector<Strand>::iterator Iter;
    for (Iter it = m_strands.begin(); it != m_strands.end(); it++)
    {
        dir = utils::genRandPointOnSphere();
        dir.m_w = 0;
        initPos = object.getMeshAABB().getCenter() + dir * object.getMeshBoundingRadius();
        dir = dir * t;

        ngl::Real length = (m_length + utils::randf(-m_lengthVariance, m_lengthVariance));
//          mass = density * volume = circle area * length = (0.25 * ngl::PI * m_thickness * m_thickness * length)
        ngl::Real strandMass = m_density * (0.25 * ngl::PI * m_thickness * m_thickness * length);
        ngl::Real pmass = strandMass / (m_nParticles - 1) + (m_nParticles - 2) * 0.5 * m_pmassFalloff;

//        initialize strand
        Strand& strand = *it;
        strand.m_sLength = length / (m_nParticles - 1);
        strand.m_pinitPos = initPos;

        initStrandParticles(strand, pmass, initPos * t, dir);
    }

    updateTransform(t);
    m_grid = new VoxelGrid(m_volume, m_volume.getWidth() / m_selfCollisionDist);
    m_grid->initialize();
}

void Hair::initStrandParticles(Strand& strand, ngl::Real pmass, ngl::Vec4 pos, const ngl::Vec4 &dir, const ngl::Vec4 &vel, const ngl::Vec4 &acc)
{
    strand.m_ppos.resize(m_nParticles);
    strand.m_pprev.resize(m_nParticles);
    strand.m_pvel.resize(m_nParticles);
    strand.m_pacc.resize(m_nParticles);
    strand.m_pmass.resize(m_nParticles);

    strand.m_ppos[0] = pos;
    strand.m_pprev[0] = pos;
    strand.m_pvel[0] = vel;
    strand.m_pacc[0] = acc;
    strand.m_pmass[0] = INFINITY;
    pos += strand.m_sLength * dir;

    for ( unsigned i = 1; i < strand.m_ppos.size(); ++i)
    {
        strand.m_ppos[i] = pos;
        strand.m_pprev[i] = pos;
        strand.m_pvel[i] = vel;
        strand.m_pacc[i] = acc;
        strand.m_pmass[i] = pmass;

        pos += strand.m_sLength * dir;
        pmass *= (1 - m_pmassFalloff);
        assert(pmass > 0);
    }
}

void Hair::updateTransform(const ngl::Mat4& t)
{
    if (m_strands.size() < 1)
        return;

    ngl::Vec4 massCenter(0,0,0,1);
    typedef std::vector<Strand>::iterator Iter;
    for (Iter it = m_strands.begin(); it != m_strands.end(); it++)
    {
        it->m_ppos[0] = it->m_pinitPos * t;
        massCenter += it->m_ppos[0];
    }
    massCenter /= m_strands.size();

    ngl::Real distSqr = 0;
    ngl::Real gridSize = 0;
    for (Iter it = m_strands.begin(); it != m_strands.end(); it++)
    {
        distSqr = (it->m_ppos[0] - massCenter).lengthSquared();
        if (distSqr > gridSize)
        {
            gridSize = distSqr;
        }
    }
    gridSize = sqrt(gridSize) + m_length;
    ngl::Vec4 vmin = massCenter - ngl::Vec4(gridSize, gridSize, gridSize, 0);
    ngl::Vec4 vmax = massCenter + ngl::Vec4(gridSize, gridSize, gridSize, 0);
    m_volume.reshape(vmin, vmax);

}

void Hair::update(ngl::Real dt)
{
    m_grid->reset();

    typedef std::vector<ngl::Vec4>::const_iterator PIter;
    typedef std::vector<Strand>::iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        Strand& strand = *it;
        updateStrand(strand, dt);
    }

//    resolveSelfInteractions(dt);
}


void Hair::resolveSelfInteractions(ngl::Real dt)
{
    typedef std::vector<Strand>::iterator SIter;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        Strand& strand = *it;
        for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
        {
            m_grid->insertDensity(strand.m_ppos[i], 1);
        }
    }

    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        Strand& strand = *it;
        for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
        {
            m_grid->insertVelocity(strand.m_ppos[i], strand.m_pvel[i]);
        }
    }

    ngl::Vec4 pressureGradient(0,0,0,0);
    ngl::Vec4 p1, p2;
    ngl::Real dr = m_selfCollisionDist;
    for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
    {
        Strand& strand = *it;
        for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
        {
//            calculate self friction
            strand.m_pvel[i] =  (1 - m_selfFriction) * strand.m_pvel[i] +
                                m_selfFriction * m_grid->getInterpolatedVelocity(strand.m_ppos[i]) / m_grid->getInterpolatedDensity(strand.m_ppos[i]);

//            calculate self repulsion
            p1.set(strand.m_ppos[i].m_x - dr, strand.m_ppos[i].m_y, strand.m_ppos[i].m_z);
            p2.set(strand.m_ppos[i].m_x + dr, strand.m_ppos[i].m_y, strand.m_ppos[i].m_z);
            pressureGradient.m_x = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

            p1.set(strand.m_ppos[i].m_x, strand.m_ppos[i].m_y - dr, strand.m_ppos[i].m_z);
            p2.set(strand.m_ppos[i].m_x, strand.m_ppos[i].m_y + dr, strand.m_ppos[i].m_z);
            pressureGradient.m_y = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

            p1.set(strand.m_ppos[i].m_x, strand.m_ppos[i].m_y, strand.m_ppos[i].m_z - dr);
            p2.set(strand.m_ppos[i].m_x, strand.m_ppos[i].m_y, strand.m_ppos[i].m_z + dr);
            pressureGradient.m_z = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

            strand.m_pvel[i] = 0.5 * (2 * strand.m_pvel[i] +  m_selfRepulsion * dt * pressureGradient / strand.m_pmass[i]);
        }
    }
}

void Hair::updateStrand(Strand &strand, ngl::Real dt)
{
    for ( unsigned i = 1; i < strand.m_ppos.size(); ++i)
    {
//        save previous position
        strand.m_pprev[i] = strand.m_ppos[i];
//        calculate new acceleration
        calculateAcceleration(strand, i, dt);
//        intergrate and calculate new position
        integrate(strand, i, dt);
    }
//    correct and update particle's positions and velocities

//    if (m_useDFTL)
//    {
////      DFTL
//        solveConstraintsDFTL(strand, dt);
//        updateVelocitiesDFTL(strand, dt);
//    } else
//    {
////      PBD
//        solveConstraints(strand, dt);
//        updateVelocities(strand, dt);
//    }


//    resolveCollisions(strand, dt);
}

void Hair::calculateAcceleration(Strand& strand, unsigned pindex, ngl::Real dt)
{
//    TODO: add other forces - wind, ect.
//    air drag force is computed by the formula 1/2 * airDensity * v^2 * areaCrossSecction * dragCoefficient
//    everything but the velocity term has been calculated and initialized in m_airDrag property in ctor
    strand.m_pacc[pindex] = utils::G - m_airDrag * strand.m_pvel[pindex].length() * strand.m_pvel[pindex] / strand.m_pmass[pindex] ;
}

void Hair::integrate(Strand &strand, unsigned pindex, ngl::Real dt)
{
//   Mid point Euler Integration
    strand.m_ppos[pindex] += dt * 0.5 * (2 * strand.m_pvel[pindex] + dt * strand.m_pacc[pindex]);
    strand.m_pvel[pindex] += dt * strand.m_pacc[pindex];
}

void Hair::solveConstraints(Strand &strand, ngl::Real dt)
{
    ngl::Real w, dist, stiffness;
    ngl::Vec4 v;
    for (unsigned k = 0; k < m_nIterations; ++k)
    {
        stiffness = 1 - pow(1 - m_stiffness, 1.0 / (k + 1));
        for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
        {
//            resolveCollisionsForParticle(strand, i);

            unsigned previ = i - 1;
            v = strand.m_ppos[previ] - strand.m_ppos[i];
            dist = v.length();
            if ( dist < utils::ERR)
                continue;

            v.normalize();
            w = (previ < 1)? 1 : (strand.m_pmass[i] / (strand.m_pmass[previ] + strand.m_pmass[i]));
            strand.m_ppos[i] += stiffness * w * (dist - strand.m_sLength) * v;
        }

        for (unsigned i = strand.m_ppos.size() - 2; i > 0; --i)
        {
//            resolveCollisionsForParticle(strand, i);

            unsigned nexti = i + 1;
            v = strand.m_ppos[nexti] - strand.m_ppos[i];
            dist = v.length();
            if ( dist < utils::ERR)
                continue;

            v.normalize();
            w = (nexti < 1)? 1 : (strand.m_pmass[i] / (strand.m_pmass[nexti] + strand.m_pmass[i]));
            strand.m_ppos[i] += stiffness * w * (dist - strand.m_sLength) * v;
        }
    }
}

void Hair::solveConstraintsDFTL(Strand &strand, ngl::Real dt)
{
    ngl::Real dist;
    ngl::Vec4 v;
    for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
    {
        resolveCollisionsForParticle(strand, i);

        unsigned previ = i - 1;
        v = strand.m_ppos[previ] - strand.m_ppos[i];
        dist = v.length();
        if ( dist < utils::ERR)
            continue;

        v.normalize();
        strand.m_pacc[i] = m_stiffness * (dist - strand.m_sLength) * v;
        strand.m_ppos[i] += strand.m_pacc[i];
    }
}

void Hair::updateVelocities(Strand &strand, ngl::Real dt)
{
    for ( unsigned i = 1; i < strand.m_ppos.size(); ++i) {
        strand.m_pvel[i] = (strand.m_ppos[i] - strand.m_pprev[i]) / dt;
    }
}

void Hair::updateVelocitiesDFTL(Strand &strand, ngl::Real dt)
{
    for ( unsigned i = 1; i < strand.m_ppos.size(); ++i) {
        ngl::Vec4 velCorrection(0,0,0,0);
        unsigned nexti = (i + 1) % strand.m_ppos.size();
        unsigned previ = std::max(i - 1, 0u);
        if (nexti > 0) {
            velCorrection += strand.m_pacc[nexti];
        } else {
            velCorrection += strand.m_pmass[previ] / (strand.m_pmass[previ] + strand.m_pmass[i]) * strand.m_pacc[i];
//            velCorrection += 0.5 * strand.m_pacc[i];
        }
        strand.m_pvel[i] = (strand.m_ppos[i] - strand.m_pprev[i] - m_damping * velCorrection)/dt;
    }
}

void Hair::resolveCollisions(Strand &strand, ngl::Real dt)
{
    std::vector<RenderObject*> obstacles;
    ngl::Vec4 collp, n, vt, vn;
    for ( unsigned i = 1; i < strand.m_ppos.size(); ++i)
    {
        m_servant->findObjectsWithinDistance( strand.m_ppos[i], strand.m_sLength + m_thickness, obstacles);
        typedef std::vector<RenderObject*>::const_iterator Iter;
        for (Iter it = obstacles.begin(); it != obstacles.end(); ++it)
        {
            if (!(*it)->isInsideObject(strand.m_ppos[i], collp, n))
                continue;

            strand.m_ppos[i] = collp + m_thickness * n;

            vn = n.dot(strand.m_pvel[i]) * n;
            vt = strand.m_pvel[i] + vn;

            if (vt.lengthSquared() > utils::ERR)
            {
                strand.m_pvel[i] -= std::max(1 - m_collisionFriction * vn.length() / vt.length(), (ngl::Real)0) * vt;
            }
        }
    }
}

void Hair::resolveCollisionsForParticle(Strand &strand, unsigned pindex)
{
    std::vector<RenderObject*> obstacles;
    ngl::Vec4 collp, n, vt, vn;

    m_servant->findObjectsWithinDistance( strand.m_ppos[pindex], strand.m_sLength + m_thickness, obstacles);
    typedef std::vector<RenderObject*>::const_iterator Iter;
    for (Iter it = obstacles.begin(); it != obstacles.end(); ++it)
    {
        if (!(*it)->isInsideObject(strand.m_ppos[pindex], collp, n))
            continue;

        strand.m_ppos[pindex] = collp + m_thickness * n;
    }
}
