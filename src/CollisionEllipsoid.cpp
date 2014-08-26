#include "CollisionEllipsoid.h"
#include "Utils.h"

CollisionEllipsoid::CollisionEllipsoid()
{
    m_localShape.identity();
    m_globalShape.identity();
    m_globalShapeInverse.identity();
}

CollisionEllipsoid::CollisionEllipsoid(const mg::Matrix4D& localShape)
{
    setLocalShape(localShape);
}

void CollisionEllipsoid::setLocalShape(const mg::Matrix4D& localShape)
{
    m_localShape = localShape;
    m_globalShape = localShape;
    m_globalShapeInverse = mg::inverse(m_globalShape);
}

void CollisionEllipsoid::updateTransform(const mg::Matrix4D& transform)
{
    m_globalShape = transform * m_localShape;
    m_globalShapeInverse = mg::inverse(m_globalShape);
}


bool CollisionEllipsoid::isInsideObject(const mg::Vec3D& p, mg::Vec3D &o_collisionPoint, mg::Vec3D &o_normal) const
{
    mg::Vec3D local_p = mg::transform_point(m_globalShapeInverse, p);
    mg::Real distance = local_p.length();
    if (distance > 1.0)
    {
        return false;
    }
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
