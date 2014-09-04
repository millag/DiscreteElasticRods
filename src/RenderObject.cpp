#include "RenderObject.h"
#include <cassert>
#include "Utils.h"

RenderObject::RenderObject(unsigned id, const Mesh *mesh):
    m_id(id), m_mesh(mesh)
{
    m_transform.identity();
    calcBoundaries();
}

RenderObject::RenderObject(const Mesh* mesh, const mg::Matrix4D &transform):
    m_id(-1), m_mesh(mesh), m_transform(transform)
{
    calcBoundaries();
}

RenderObject::RenderObject(unsigned id, const Mesh* mesh, const mg::Matrix4D& transform):
    m_id(id), m_mesh(mesh), m_transform(transform)
{
    calcBoundaries();
}




void RenderObject::setTransform(const mg::Matrix4D &t)
{
    m_transform = t;
    typedef std::vector<CollisionShape>::iterator Iter;
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
    m_meshBoundingRadius = m_meshAABB.getEnclosedRadius();

    m_boundingRadius = mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EX).length();
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EY).length());
    m_boundingRadius = std::max(m_boundingRadius, mg::transform_vector(m_transform, m_meshBoundingRadius * mg::EZ).length());
}

void RenderObject::addCollisionShape(const CollisionShape& shape)
{
    unsigned idx = m_collisionShapes.size();
    m_collisionShapes.push_back(shape);
    m_collisionShapes[idx].updateTransform(m_transform);
}

bool RenderObject::isInsideObject(const mg::Vec3D& p, mg::Vec3D &o_collisionPoint, mg::Vec3D &o_normal) const
{
    typedef std::vector<CollisionShape>::const_iterator Iter;
    for (Iter it = m_collisionShapes.begin(); it != m_collisionShapes.end(); ++it)
    {
        if (it->isInside(p, o_collisionPoint, o_normal))
        {
            return true;
        }
    }
    return false;
}
