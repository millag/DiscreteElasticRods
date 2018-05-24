#ifndef STRAND_H
#define STRAND_H

#include <vector>
#include "IServant.h"
#include "VoxelGrid.h"

// length - measured in meters 1gl unit = 1 m
// mass - measured in kg
// time - measured in sec
// velocity - m/s
// acceleration - m/s^2
// density - mass per volume kg / m^3



struct Strand
{
public:
    ngl::Real m_sLength;
    ngl::Vec4 m_pinitPos;
    std::vector<ngl::Vec4> m_ppos;
    std::vector<ngl::Vec4> m_pprev;
    std::vector<ngl::Vec4> m_pvel;
    std::vector<ngl::Vec4> m_pacc;
    std::vector<ngl::Real> m_pmass;
};


class Hair
{
public:
    Hair(IServant* servant);
    ~Hair();

    void initialize(const RenderObject &object, unsigned nMaxStrands);
    void updateTransform(const ngl::Mat4& t);
    void update(ngl::Real dt);

    std::vector<Strand> m_strands;

    ngl::Real m_length;
    ngl::Real m_lengthVariance;
    ngl::Real m_density;
    ngl::Real m_pmassFalloff;
    ngl::Real m_thickness;
    ngl::Real m_stiffness;
    ngl::Real m_damping;
    ngl::Real m_collisionFriction;
    ngl::Real m_airDrag;
    ngl::Real m_selfCollisionDist;
    ngl::Real m_selfFriction;
    ngl::Real m_selfRepulsion;
    bool m_useDFTL;

    unsigned m_nParticles;
    unsigned m_nIterations;

private:
    Hair():m_servant(NULL) { }

    IServant* m_servant;
    VoxelGrid* m_grid;
    AABB m_volume;

    void initStrandParticles(Strand& strand, ngl::Real pmass, ngl::Vec4 pos,
                   const ngl::Vec4 &dir = ngl::Vec4(0,-1,0,0),
                   const ngl::Vec4 &vel = ngl::Vec4(0,0,0,0), const ngl::Vec4 &acc = ngl::Vec4(0,0,0,0));

    void updateStrand(Strand &strand, ngl::Real dt);
    void calculateAcceleration(Strand &strand, unsigned pindex, ngl::Real dt);
    void integrate(Strand &strand, unsigned pindex, ngl::Real dt);
    void solveConstraintsDFTL(Strand &strand, ngl::Real dt);
    void updateVelocitiesDFTL(Strand &strand, ngl::Real dt);
    void solveConstraints(Strand &strand, ngl::Real dt);
    void updateVelocities(Strand &strand, ngl::Real dt);
    void resolveCollisions(Strand &strand, ngl::Real dt);
    void resolveCollisionsForParticle(Strand &strand, unsigned pindex);
    void resolveSelfInteractions(ngl::Real dt);
    void updateVolume();

    ngl::Vec4 getAirDrag(Strand &strand, unsigned pindex);
};

#endif // STRAND_H
