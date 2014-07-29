#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <ngl/Mat4.h>
#include "Mesh.h"
#include "AABB.h"


class RenderObject
{

public:
    RenderObject(const Mesh* _mesh = NULL, const ngl::Mat4& _transform = ngl::Mat4(), int _shaderId = -1);
    virtual ~RenderObject() {}

    unsigned getMeshId() const;
    const Mesh* getMesh() const;
    void setMesh(const Mesh* _mesh);

    virtual const ngl::Mat4& getTransform() const { return m_transform; }
    virtual void setTransform(const ngl::Mat4& t);
    ngl::Vec4 getPosition() const  { return m_AABB.getCenter(); }

    ngl::Real getBoundingRadius() const { return m_boundingRadius; }
    const AABB& getAABB() const { return m_AABB; }

    ngl::Real getMeshBoundingRadius() const { return m_meshAABB.getBoundingRadius(); }
    const AABB& getMeshAABB() const { return m_meshAABB; }

    bool isInsideObject(const ngl::Vec4& p, ngl::Vec4& o_p, ngl::Vec4& o_n) const;

protected:
    void calcBoundaries();
    void calcAABB();

    const Mesh* m_mesh;
    ngl::Mat4 m_transform;
    ngl::Real m_boundingRadius;
    AABB m_meshAABB;
    AABB m_AABB;

    int m_shaderId;
};

#endif // GEOMETRY_H
