#pragma once

#include "ElasticRodState.h"
#include <vector>

class Hair;

struct HairState
{
public:
	inline void clear()
	{
		m_strands.clear();
	}

private:
	std::vector<ElasticRodState> m_strands;

private:
	friend class Hair;
};
