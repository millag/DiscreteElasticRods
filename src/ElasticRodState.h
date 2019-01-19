#pragma once

#include "config.h"

class ElasticRod;

struct ElasticRodState
{
public:
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
