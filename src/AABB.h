#pragma once

#include "config.h"

class AABB
{
public:
	AABB():
	    m_vmin(0,0,0)
	  , m_vmax(0,0,0)
	{}

	AABB(const mg::Vec3D& p1,  const mg::Vec3D& p2):
	    m_vmin(p1)
	  , m_vmax(p2)
	{
		adjustBounds();
	}

	inline bool isEmpty() const
	{
		return m_vmax[0] <= m_vmin[0] || m_vmax[1] <= m_vmin[1] || m_vmax[2] <= m_vmin[2];
	}

	inline void setEmpty()
	{
		m_vmax = m_vmin = mg::Origin;
	}

	inline const mg::Vec3D& getMin() const { return m_vmin; }
	inline const mg::Vec3D& getMax() const { return m_vmax; }
	inline mg::Vec3D getCenter() const { return (m_vmin + m_vmax) * 0.5; }
	inline mg::Vec3D getSize() const { return m_vmax - m_vmin; }
	inline mg::Real getWidth() const { return getSize()[0]; }
	inline mg::Real getHeight() const { return getSize()[1]; }
	inline mg::Real getDepth() const { return getSize()[2]; }

	inline mg::Real getBoundingRadius() const
	{
		return getSize().length() * 0.5f;
	}

	inline mg::Real getEnclosedRadius() const
	{
		auto d = getSize();
		return std::min(d[0], std::min(d[1], d[2])) * 0.5f;
	}

	inline mg::Vec3D getBLF() const { return m_vmin; }
	inline mg::Vec3D getBLB() const { return mg::Vec3D(m_vmin[0], m_vmin[1], m_vmax[2]); }
	inline mg::Vec3D getBRF() const { return mg::Vec3D(m_vmax[0], m_vmin[1], m_vmin[2]); }
	inline mg::Vec3D getBRB() const { return mg::Vec3D(m_vmax[0], m_vmin[1], m_vmax[2]); }

	inline mg::Vec3D getTRB() const { return m_vmax; }
	inline mg::Vec3D getTRF() const { return mg::Vec3D(m_vmax[0], m_vmax[1], m_vmin[2]); }
	inline mg::Vec3D getTLB() const { return mg::Vec3D(m_vmin[0], m_vmax[1], m_vmax[2]); }
	inline mg::Vec3D getTLF() const { return mg::Vec3D(m_vmin[0], m_vmax[1], m_vmin[2]); }

	inline void reshape(const mg::Vec3D& p1, const mg::Vec3D& p2)
	{
		m_vmin = p1;
		m_vmax = p2;
		adjustBounds();
	}

protected:
	inline void adjustBounds()
	{
		for (auto i = 0; i < 3; ++i)
		{
			if (m_vmax[i] < m_vmin[i])
			{
				std::swap(m_vmax[i], m_vmin[i]);
			}
		}
	}

	mg::Vec3D m_vmin;
	mg::Vec3D m_vmax;
};
