#ifndef HAIRSTRAND_H
#define HAIRSTRAND_H

#include "ElasticRod.h"
#include "RenderObject.h"

class HairStrand
{
public:
    HairStrand();
    ~HairStrand();

    const std::vector<ElasticRod*>& getStrands() const { return m_strands; }

    void initialize(const RenderObject* object);
    void update(mg::Real dt);

    mg::Real m_length;
    mg::Real m_lengthVariance;
    mg::Real m_density;
    mg::Real m_thickness;
    mg::Real m_bendStiffness;
    mg::Real m_twistStiffness;
    mg::Real m_maxForce;

    unsigned m_nParticles;
    unsigned m_nIterations;

private:
    std::vector<ElasticRod*> m_strands;
    const RenderObject* m_object;

};

#endif // HAIRSTRAND_H
