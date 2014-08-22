#ifndef HAIR_H
#define HAIR_H

#include "ElasticRod.h"
#include "RenderObject.h"

struct HairParams
{
    HairParams();
    ~HairParams();

    mg::Real m_length;
    mg::Real m_lengthVariance;
    mg::Real m_helicalRadius;
    mg::Real m_helicalPitch;

    mg::Real m_density;
    mg::Real m_thickness;
    unsigned m_nParticles;

    RodParams* m_rodParams;
};

class Hair
{
public:
    Hair();
    ~Hair();

    void reset();
    void update(mg::Real dt);

public:
    HairParams* m_params;
    std::vector<ElasticRod*> m_strands;
    const RenderObject* m_object;
    std::vector<unsigned> m_findices;
    std::vector<unsigned> m_vindices;
};

#endif // HAIR_H
