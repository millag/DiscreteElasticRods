#include "AABB.h"
#include "Utils.h"

AABB::AABB():m_vmin(),m_vmax(),m_boundingRadius(0)
{ }

AABB::AABB(const ngl::Vec4 &_vmin, const ngl::Vec4 &_vmax):m_vmin(_vmin),m_vmax(_vmax),m_boundingRadius(0)
{
    assert(m_vmin.m_x <= m_vmax.m_x && m_vmin.m_y <= m_vmax.m_y && m_vmin.m_z <= m_vmax.m_z);
    m_boxSize = m_vmax - m_vmin;
    calcBoundaries();
}

void AABB::reshape(const ngl::Vec4 &_vmin, const ngl::Vec4 &_vmax)
{
    assert(_vmin.m_x <= _vmax.m_x && _vmin.m_y <= _vmax.m_y && _vmin.m_z <= _vmax.m_z);
    m_vmin = _vmin;
    m_vmax = _vmax;
    calcBoundaries();
}

void AABB::calcBoundaries()
{
    m_boxSize = m_vmax - m_vmin;
   // m_boundingRadius = (m_vmax - getCenter()).length();
    m_boundingRadius = std::max(m_boxSize.m_x, std::max(m_boxSize.m_y, m_boxSize.m_z)) / 2;
    assert(fabs((m_vmax - getCenter()).length() - (m_vmin - getCenter()).length()) < mg::ERR);
}
