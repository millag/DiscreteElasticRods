#ifndef AABB_H
#define AABB_H

#include "config.h"

class AABB
{
public:
    AABB();
    AABB(const mg::Vec3D& _vmin,  const mg::Vec3D& _vmax);

    inline mg::Vec3D getVMin() const { return m_vmin; }
    inline mg::Vec3D getVMax() const { return m_vmax; }
    inline mg::Vec3D getCenter() const { return (m_vmin + m_vmax) * 0.5; }

    inline mg::Vec3D getBLF() const { return m_vmin; }
    inline mg::Vec3D getBLB() const { return mg::Vec3D(m_vmin[0], m_vmin[1], m_vmax[2]); }
    inline mg::Vec3D getBRF() const { return mg::Vec3D(m_vmax[0], m_vmin[1], m_vmin[2]); }
    inline mg::Vec3D getBRB() const { return mg::Vec3D(m_vmax[0], m_vmin[1], m_vmax[2]); }

    inline mg::Vec3D getTRB() const { return m_vmax; }
    inline mg::Vec3D getTRF() const { return mg::Vec3D(m_vmax[0], m_vmax[1], m_vmin[2]); }
    inline mg::Vec3D getTLB() const { return mg::Vec3D(m_vmin[0], m_vmax[1], m_vmax[2]); }
    inline mg::Vec3D getTLF() const { return mg::Vec3D(m_vmin[0], m_vmax[1], m_vmin[2]); }

    inline mg::Real getWidth() const { return m_boxSize[0]; }
    inline mg::Real getHeight() const { return m_boxSize[1]; }
    inline mg::Real getDepth() const { return m_boxSize[2]; }
    inline mg::Real getEnclosedRadius() const { return m_enclosedRadius; }

    void reshape(const mg::Vec3D& _vmin,  const mg::Vec3D& _vmax);

protected:
    void calcBoundaries();

    mg::Vec3D m_vmin;
    mg::Vec3D m_vmax;
    mg::Vec3D m_boxSize;
    mg::Real m_enclosedRadius;

};


#endif // AABB_H
