#ifndef SPIRAL_H
#define SPIRAL_H

#include "RenderObject.h"
#include "ElasticRod.h"

class Spiral
{
public:
    Spiral();
    ~Spiral();

    const std::vector<ElasticRod*>& getStrands() const { return m_strands; }

    void init(const RenderObject* object);
    void update(mg::Real dt);

    mg::Real m_radius;
    mg::Real m_lenght;

    mg::Real m_bendStiffness;
    mg::Real m_twistStiffness;
    mg::Real m_maxForce;

    unsigned m_nParticles;
    unsigned m_nIterations;

private:
    std::vector<ElasticRod*> m_strands;
    const RenderObject* m_object;
};

#endif // SPIRAL_H
