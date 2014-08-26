#ifndef COLLISIONOBJECT_H
#define COLLISIONOBJECT_H

#include "config.h"

class CollisionEllipsoid
{
public:
    CollisionEllipsoid();
    CollisionEllipsoid(const mg::Matrix4D& localShape);
    ~CollisionEllipsoid() { }

    void setLocalShape(const mg::Matrix4D& localShape);
    void updateTransform(const mg::Matrix4D& transform);

    bool isInsideObject(const mg::Vec3D& p, mg::Vec3D& o_collisionPoint, mg::Vec3D& o_normal) const;

protected:
    mg::Matrix4D m_localShape;
    mg::Matrix4D m_globalShape;
    mg::Matrix4D m_globalShapeInverse;
};
#endif // COLLISIONOBJECT_H
