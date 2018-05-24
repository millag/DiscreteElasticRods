#ifndef HAIRSTATE_H
#define HAIRSTATE_H

#include <vector>
#include "ElasticRodState.h"

class Hair;

struct HairState
{
public:
    HairState() { clear(); }
    ~HairState() { clear(); }

    inline void clear() { m_strands.clear(); }

private:
    std::vector<ElasticRodState> m_strands;

private:
    friend class Hair;
};

#endif // HAIRSTATE_H
