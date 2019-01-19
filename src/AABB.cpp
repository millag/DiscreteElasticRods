#include "AABB.h"
#include <cassert>
#include "Utils.h"

AABB::AABB(const mg::Vec3D& p1, const mg::Vec3D& p2):
    m_vmin(p1)
  , m_vmax(p2)
{
	void adjustBounds();
}

void AABB::setEmpty()
{
	m_vmax = m_vmin = mg::Origin;
}

void AABB::reshape(const mg::Vec3D& p1, const mg::Vec3D& p2)
{
	m_vmin = p1;
	m_vmax = p2;
	adjustBounds();
}

void AABB::adjustBounds()
{
	for (auto i = 0; i < 3; ++i)
	{
		if (m_vmax[i] < m_vmin[i])
		{
			std::swap(m_vmax[i], m_vmin[i]);
		}
	}
}
