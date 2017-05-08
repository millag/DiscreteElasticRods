#include "AABB.h"
#include <cassert>
#include "Utils.h"


AABB::AABB():m_vmin(0,0,0),m_vmax(0,0,0),m_enclosedRadius(0)
{ }

AABB::AABB(const mg::Vec3D &_vmin, const mg::Vec3D &_vmax):m_vmin(_vmin),m_vmax(_vmax),m_enclosedRadius(0)
{
    assert(_vmin[0] <= _vmax[0] && _vmin[1] <= _vmax[1] && _vmin[2] <= _vmax[2]);
    m_boxSize = m_vmax - m_vmin;
    calcBoundaries();
}

void AABB::reshape(const mg::Vec3D &_vmin, const mg::Vec3D &_vmax)
{
    assert(_vmin[0] <= _vmax[0] && _vmin[1] <= _vmax[1] && _vmin[2] <= _vmax[2]);
    m_vmin = _vmin;
    m_vmax = _vmax;
    calcBoundaries();
}

void AABB::calcBoundaries()
{
    m_boxSize = m_vmax - m_vmin;
   // m_boundingRadius = (m_vmax - getCenter()).length();
    m_enclosedRadius = std::max(m_boxSize[0], std::max(m_boxSize[1], m_boxSize[2])) * 0.5f;
    assert(fabs((m_vmax - getCenter()).length() - (m_vmin - getCenter()).length()) < mg::ERR);
}
