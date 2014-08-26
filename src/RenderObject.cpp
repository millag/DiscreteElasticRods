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
    typedef std::vector<CollisionEllipsoid>::iterator Iter;
    for (Iter it = m_collisionShapes.begin(); it != m_collisionShapes.end(); ++it)
    {
        it->updateTransform(m_transform);
    }

    m_boundingRadius = mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EX).length();
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EY).length());
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EZ).length());
}

void RenderObject::calcBoundaries()
{

    if (!m_mesh || !m_mesh->m_vertices.size())
    {
        m_meshAABB.reshape(mg::Vec3D(0,0,0), mg::Vec3D(0,0,0));
        return;
    }

    mg::Vec3D vmin = *(m_mesh->m_vertices.begin());
    mg::Vec3D vmax = vmin;
    m_meshBoundingRadius = vmin.length_squared();

    typedef std::vector<mg::Vec3D>::const_iterator VIter;
    for (VIter it = m_mesh->m_vertices.begin() + 1; it != m_mesh->m_vertices.end(); ++it) {
        const mg::Vec3D &vert = (*it);
        vmin[0] = std::min(vert[0], vmin[0]);
        vmin[1] = std::min(vert[1], vmin[1]);
        vmin[2] = std::min(vert[2], vmin[2]);

        vmax[0] = std::max(vert[0], vmax[0]);
        vmax[1] = std::max(vert[1], vmax[1]);
        vmax[2] = std::max(vert[2], vmax[2]);

        mg::Real lengthSqr = vert.length_squared();
        if (m_meshBoundingRadius < lengthSqr)
        {
            m_meshBoundingRadius = lengthSqr;
        }
    }

    m_meshAABB.reshape(vmin, vmax);
    m_meshBoundingRadius = std::sqrt(m_meshBoundingRadius);

    m_boundingRadius = mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EX).length();
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EY).length());
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EZ).length());
}

void RenderObject::addCollisionShape(const CollisionEllipsoid& ellipsoid)
{
    unsigned idx = m_collisionShapes.size();
    m_collisionShapes.push_back(ellipsoid);
    m_collisionShapes[idx].updateTransform(m_transform);
}

bool RenderObject::isInsideObject(const mg::Vec3D& p, mg::Vec3D &o_collisionPoint, mg::Vec3D &o_normal) const
{
    typedef std::vector<CollisionEllipsoid>::const_iterator Iter;
    for (Iter it = m_collisionShapes.begin(); it != m_collisionShapes.end(); ++it)
    {
        if (it->isInsideObject(p, o_collisionPoint, o_normal))
        {
            return true;
        }
    }
    return false;
}
