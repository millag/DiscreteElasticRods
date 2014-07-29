#ifndef AABB_H
#define AABB_H

#include <ngl/Vec4.h>

class AABB
{
public:
    AABB();
    AABB(const ngl::Vec4& _vmin,  const ngl::Vec4& _vmax);

    ngl::Vec4 getVMin() const { return m_vmin; }
    ngl::Vec4 getVMax() const { return m_vmax; }
    ngl::Vec4 getCenter() const { return (m_vmin + m_vmax) / 2; }

    ngl::Vec4 getBLF() const { return m_vmin; }
    ngl::Vec4 getBLB() const { return ngl::Vec4(m_vmin.m_x, m_vmin.m_y, m_vmax.m_z); }
    ngl::Vec4 getBRF() const { return ngl::Vec4(m_vmax.m_x, m_vmin.m_y, m_vmin.m_z); }
    ngl::Vec4 getBRB() const { return ngl::Vec4(m_vmax.m_x, m_vmin.m_y, m_vmax.m_z); }

    ngl::Vec4 getTRB() const { return m_vmax; }
    ngl::Vec4 getTRF() const { return ngl::Vec4(m_vmax.m_x, m_vmax.m_y, m_vmin.m_z); }
    ngl::Vec4 getTLB() const { return ngl::Vec4(m_vmin.m_x, m_vmax.m_y, m_vmax.m_z); }
    ngl::Vec4 getTLF() const { return ngl::Vec4(m_vmin.m_x, m_vmax.m_y, m_vmin.m_z); }

    ngl::Real getWidth() const { return m_boxSize.m_x; }
    ngl::Real getHeight() const { return m_boxSize.m_y; }
    ngl::Real getDepth() const { return m_boxSize.m_z; }
    ngl::Real getBoundingRadius() const { return m_boundingRadius; }

    void reshape(const ngl::Vec4& _vmin,  const ngl::Vec4& _vmax);

protected:
    void calcBoundaries();

    ngl::Vec4 m_vmin;
    ngl::Vec4 m_vmax;
    ngl::Vec4 m_boxSize;
    ngl::Real m_boundingRadius;

};


#endif // AABB_H
