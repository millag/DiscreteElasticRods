#include "RenderObject.h"
#include <cassert>
#include "Utils.h"

RenderObject::RenderObject(const Mesh* _mesh, const mg::Matrix4D &_transform, int _shaderId):
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

void RenderObject::setTransform(const mg::Matrix4D &t)
{
    m_transform = t;
    calcAABB();
}

void RenderObject::calcBoundaries()
{

    if (!m_mesh || !m_mesh->m_vertices.size())
    {
        m_meshAABB.reshape(mg::Vec3D(), mg::Vec3D());
        m_AABB.reshape(mg::Vec3D(), mg::Vec3D());
        return;
    }

    mg::Vec3D vmin = *(m_mesh->m_vertices.begin());
    mg::Vec3D vmax = vmin;
    typedef std::vector<mg::Vec3D>::const_iterator VIter;
    for (VIter it = m_mesh->m_vertices.begin() + 1; it != m_mesh->m_vertices.end(); ++it) {
        const mg::Vec3D &vert = (*it);
        vmin[0] = std::min(vert[0], vmin[0]);
        vmin[1] = std::min(vert[1], vmin[1]);
        vmin[2] = std::min(vert[2], vmin[2]);

        vmax[0] = std::max(vert[0], vmax[0]);
        vmax[1] = std::max(vert[1], vmax[1]);
        vmax[2] = std::max(vert[2], vmax[2]);
    }

    m_meshAABB.reshape(vmin, vmax);
//  apply current transform and recalculate AABB
    calcAABB();
    m_boundingRadius = m_AABB.getBoundingRadius();
}

//    TODO FIX need to add transformation somehow
void RenderObject::calcAABB()
{

    std::vector<mg::Vec3D> corners;
    corners.reserve(8);
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getBLF()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getBLB()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getBRF()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getBRB()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getTRB()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getTRF()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getTLB()) );
    corners.push_back( mg::transform_point(m_transform, m_meshAABB.getTLF()) );

    mg::Vec3D vmin = *corners.begin();
    mg::Vec3D vmax = vmin;
    typedef std::vector<mg::Vec3D>::iterator VIter;
    for (VIter it = corners.begin() + 1; it != corners.end(); ++it) {
        const mg::Vec3D &vert = (*it);
        vmin[0] = std::min(vert[0], vmin[0]);
        vmin[1] = std::min(vert[1], vmin[1]);
        vmin[2] = std::min(vert[2], vmin[2]);

        vmax[0] = std::max(vert[0], vmax[0]);
        vmax[1] = std::max(vert[1], vmax[1]);
        vmax[2] = std::max(vert[2], vmax[2]);
    }

    m_AABB.reshape(vmin, vmax);
}

bool RenderObject::isInsideObject(const mg::Vec3D& p, mg::Vec3D &o_p, mg::Vec3D &o_n) const
{
    o_n = p - getPosition();
    mg::Real distSqr = o_n.length_squared();

    if (distSqr > getBoundingRadius() * getBoundingRadius())
    {
        return false;
    }
    if (distSqr > mg::ERR)
    {
        o_n.normalize();
    }

    o_p =  getPosition() + o_n  * getBoundingRadius();
    return true;
}
