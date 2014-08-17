#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Mesh.h"
#include "AABB.h"


class RenderObject
{

public:
    RenderObject(const Mesh* _mesh, const mg::Matrix4D& _transform , int _shaderId = -1);
    virtual ~RenderObject() {}

    unsigned getMeshId() const;
    const Mesh* getMesh() const;
    void setMesh(const Mesh* _mesh);

    virtual const mg::Matrix4D& getTransform() const { return m_transform; }
    virtual void setTransform(const mg::Matrix4D& t);
    mg::Vec3D getPosition() const  { return m_AABB.getCenter(); }

    mg::Real getBoundingRadius() const { return m_boundingRadius; }
    const AABB& getAABB() const { return m_AABB; }

    mg::Real getMeshBoundingRadius() const { return m_meshAABB.getBoundingRadius(); }
    const AABB& getMeshAABB() const { return m_meshAABB; }

    bool isInsideObject(const mg::Vec3D& p, mg::Vec3D& o_p, mg::Vec3D& o_n) const;

protected:
    void calcBoundaries();
    void calcAABB();

    const Mesh* m_mesh;
    mg::Matrix4D m_transform;
    mg::Real m_boundingRadius;
    AABB m_meshAABB;
    AABB m_AABB;

    int m_shaderId;
};

#endif // GEOMETRY_H
