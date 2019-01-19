#pragma once

#include "Hair.h"
#include <deque>

class AnimationBuffer
{
public:

	AnimationBuffer(unsigned capacity = 1000):
	    m_capacity(capacity)
	{ }

	inline void setCapacity(unsigned capacity)
	{
		m_capacity = capacity;
		m_buffer.clear();
		m_state.clear();
	}

	inline unsigned getCapacity() const  { return m_capacity; }
	inline void clear() { m_state.clear(); m_buffer.clear(); }
	inline unsigned size() const { return m_buffer.size(); }
	inline const mg::Matrix4D& getFrame(unsigned idx) const { return m_buffer[idx]; }

	inline void saveHairState(const Hair* hair) { hair->getState(m_state); }
	inline void restoreHairState(Hair* hair) { hair->setState(m_state); }

	inline void saveFrame(const mg::Matrix4D transform)
	{
		if (m_buffer.size() == m_capacity)
		{
			return;
		}
		m_buffer.push_back(transform);
	}


private:

	unsigned  m_capacity;
	std::deque<mg::Matrix4D> m_buffer;
	HairState m_state;
};
