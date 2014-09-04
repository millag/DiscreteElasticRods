#ifndef ELASTICRODSTATE_H
#define ELASTICRODSTATE_H

#include "config.h"

class ElasticRod;

struct ElasticRodState
{
public:
    ElasticRodState() { clear(); }
    ~ElasticRodState() {  clear(); }

    inline void clear()
    {
        m_ppos.clear();
        m_pvel.clear();
        m_u0.zero();
    }

private:
    std::vector<mg::Vec3D> m_ppos;
    std::vector<mg::Vec3D> m_pvel;
    mg::Vec3D m_u0;

private:
    friend class ElasticRod;
};

#endif // ELASTICRODSTATE_H
