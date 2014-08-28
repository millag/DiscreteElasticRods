#include "CollisionShape.h"
#include "Utils.h"

CollisionShape::CollisionShape()
{
    m_localShape.identity();
    m_globalShape.identity();
    m_globalShapeInverse.identity();
}

CollisionShape::CollisionShape(const mg::Matrix4D& localShape)
{
    setLocalShape(localShape);
}

void CollisionShape::setLocalShape(const mg::Matrix4D& localShape)
{
    m_localShape = localShape;
    m_globalShape = localShape;
    m_globalShapeInverse = mg::inverse(m_globalShape);
}

void CollisionShape::updateTransform(const mg::Matrix4D& transform)
{
    m_globalShape = transform * m_localShape;
    m_globalShapeInverse = mg::inverse(m_globalShape);
}

// Citation Begin
// code borrowed and modified from http://code.google.com/p/opencloth/

bool CollisionShape::isInside(const mg::Vec3D& p, mg::Vec3D &o_collisionPoint, mg::Vec3D &o_normal) const
{
    mg::Vec3D local_p = mg::transform_point(m_globalShapeInverse, p);
    mg::Real distance = local_p.length_squared();
    if (distance > 1.0)
    {
        return false;
    }
    distance = std::sqrt(distance);
    mg::Vec3D local_dir = (1.0 - distance) * local_p / distance;

    // Transform back to original space
    mg::Vec3D transformInv;
    transformInv.set(m_globalShape(0, 0), m_globalShape(0, 1), m_globalShape(0, 2));
    o_normal[0] = mg::dot(transformInv / mg::dot(transformInv, transformInv), local_dir);

    transformInv.set(m_globalShape(1, 0), m_globalShape(1, 1), m_globalShape(1, 2));
    o_normal[1] = mg::dot(transformInv / mg::dot(transformInv, transformInv), local_dir);

    transformInv.set(m_globalShape(2, 0), m_globalShape(2, 1), m_globalShape(2, 2));
    o_normal[2] = mg::dot(transformInv / mg::dot(transformInv, transformInv), local_dir);
    o_collisionPoint = p + o_normal;

    return true;
}

// Citation End
