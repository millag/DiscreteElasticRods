#ifndef COLLISIONSHAPE_H
#define COLLISIONSHAPE_H

#include "config.h"

class CollisionShape
{
public:
    CollisionShape();
    CollisionShape(const mg::Matrix4D& localShape);
    ~CollisionShape() { }

    void setLocalShape(const mg::Matrix4D& localShape);
    void updateTransform(const mg::Matrix4D& transform);

    bool isInside(const mg::Vec3D& p, mg::Vec3D& o_collisionPoint, mg::Vec3D& o_normal) const;

protected:
    mg::Matrix4D m_localShape;
    mg::Matrix4D m_globalShape;
    mg::Matrix4D m_globalShapeInverse;
};
#endif // COLLISIONSHAPE_H
