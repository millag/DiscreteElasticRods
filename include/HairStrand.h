#ifndef HAIRSTRAND_H
#define HAIRSTRAND_H

#include "ElasticRod.h"
#include "RenderObject.h"


struct HairParams
{
    mg::Real m_length;
    mg::Real m_lengthVariance;
    mg::Real m_density;
    mg::Real m_thickness;
    mg::Real m_bendStiffness;
    mg::Real m_twistStiffness;
    mg::Real m_maxForce;

    unsigned m_nParticles;
    unsigned m_nIterations;
};

class HairStrand
{
public:
    HairStrand();
    ~HairStrand();

    const std::vector<ElasticRod*>& getStrands() const { return m_strands; }

    void generateStrands(const RenderObject* object);
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
    RodParams* m_rodParams;
    const RenderObject* m_object;

};

#endif // HAIRSTRAND_H
