#include "RenderObject.h"
#include "Utils.h"

RenderObject::RenderObject(const Mesh* _mesh, const ngl::Mat4 &_transform, int _shaderId):
    m_mesh(_mesh), m_transform(_transform), m_shaderId(_shaderId)
{
    calcBoundaries();
}


unsigned RenderObject::getMeshId() const
{
    assert(m_mesh != NULL);
    return m_mesh->getId();
}

const Mesh* RenderObject::getMesh() const
{
    return m_mesh;
}

void RenderObject::setMesh(const Mesh* _mesh)
{
    m_mesh = _mesh;
    calcBoundaries();
}

void RenderObject::setTransform(const ngl::Mat4 &t)
{
    m_transform = t;
    calcAABB();
}

void RenderObject::calcBoundaries()
{

    if (!m_mesh || !m_mesh->m_vertices.size()) {
        m_meshAABB.reshape(ngl::Vec4(), ngl::Vec4());
        m_AABB.reshape(ngl::Vec4(), ngl::Vec4());
        return;
    }

    ngl::Vec4 vmin = *(m_mesh->m_vertices.begin());
    ngl::Vec4 vmax = vmin;
    typedef std::vector<ngl::Vec4>::const_iterator VIter;
    for (VIter it = m_mesh->m_vertices.begin() + 1; it != m_mesh->m_vertices.end(); ++it) {
        const ngl::Vec4 &vert = (*it);
        vmin.m_x = std::min(vert.m_x, vmin.m_x);
        vmin.m_y = std::min(vert.m_y, vmin.m_y);
        vmin.m_z = std::min(vert.m_z, vmin.m_z);

        vmax.m_x = std::max(vert.m_x, vmax.m_x);
        vmax.m_y = std::max(vert.m_y, vmax.m_y);
        vmax.m_z = std::max(vert.m_z, vmax.m_z);
    }

    m_meshAABB.reshape(vmin, vmax);
//  apply current transform and recalculate AABB
    calcAABB();
    m_boundingRadius = m_AABB.getBoundingRadius();
}

void RenderObject::calcAABB()
{
    std::vector<ngl::Vec4> corners;
    corners.reserve(8);
    corners.push_back(m_meshAABB.getBLF() * m_transform);
    corners.push_back(m_meshAABB.getBLB() * m_transform);
    corners.push_back(m_meshAABB.getBRF() * m_transform);
    corners.push_back(m_meshAABB.getBRB() * m_transform);
    corners.push_back(m_meshAABB.getTRB() * m_transform);
    corners.push_back(m_meshAABB.getTRF() * m_transform);
    corners.push_back(m_meshAABB.getTLB() * m_transform);
    corners.push_back(m_meshAABB.getTLF() * m_transform);

    ngl::Vec4 vmin = *corners.begin();
    ngl::Vec4 vmax = vmin;
    typedef std::vector<ngl::Vec4>::iterator VIter;
    for (VIter it = corners.begin() + 1; it != corners.end(); ++it) {
        const ngl::Vec4 &vert = (*it);
        vmin.m_x = std::min(vert.m_x, vmin.m_x);
        vmin.m_y = std::min(vert.m_y, vmin.m_y);
        vmin.m_z = std::min(vert.m_z, vmin.m_z);

        vmax.m_x = std::max(vert.m_x, vmax.m_x);
        vmax.m_y = std::max(vert.m_y, vmax.m_y);
        vmax.m_z = std::max(vert.m_z, vmax.m_z);
    }

    m_AABB.reshape(vmin, vmax);
}

bool RenderObject::isInsideObject(const ngl::Vec4& p, ngl::Vec4 &o_p, ngl::Vec4 &o_n) const
{
    o_n = p - getPosition();
    ngl::Real distSqr = o_n.lengthSquared();

    if (distSqr > getBoundingRadius() * getBoundingRadius())
    {
        return false;
    }

    if (distSqr < utils::ERR)
    {
        o_n = utils::genRandPointOnSphere() ;
        o_n.m_w = 0;
    }
    o_n.normalize();
    o_p =  getPosition() + o_n  * getBoundingRadius();
    return true;
}
