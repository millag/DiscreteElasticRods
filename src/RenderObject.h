#ifndef RENDEROBJECT_H
#define RENDEROBJECT_H

#include "Mesh.h"
#include "AABB.h"
#include "CollisionShape.h"


class RenderObject
{
public:

    RenderObject(unsigned id = -1, const Mesh* mesh = NULL);
    RenderObject(const Mesh* mesh, const mg::Matrix4D& transform);
    RenderObject(unsigned id, const Mesh* mesh, const mg::Matrix4D& transform);

    ~RenderObject() {}

    inline unsigned getId() const { return m_id; }
    inline void setId(unsigned id) { m_id = id; }

    inline unsigned getMeshId() const { assert(m_mesh != NULL); return m_mesh->getId();}
    inline const Mesh* getMesh() const { return m_mesh; }
    inline void setMesh(const Mesh* _mesh) { m_mesh = _mesh; calcBoundaries();}

    void setTransform(const mg::Matrix4D& t);
    inline const mg::Matrix4D& getTransform() const { return m_transform; }
    inline mg::Vec3D getPosition() const  { return mg::matrix_get_translation(m_transform); }

    inline mg::Vec3D getCenter() const  { return mg::transform_point(m_transform, m_meshAABB.getCenter()); }
    inline mg::Real getBoundingRadius() const { return m_boundingRadius; }
    inline mg::Real getMeshBoundingRadius() const { return m_meshBoundingRadius; }
    inline const AABB& getMeshAABB() const { return m_meshAABB; }

    void addCollisionShape(const CollisionShape& shape);
    bool isInsideObject(const mg::Vec3D& p, mg::Vec3D& o_collisionPoint, mg::Vec3D& o_normal) const;

protected:

    void calcBoundaries();

    unsigned m_id;

    const Mesh* m_mesh;
    mg::Matrix4D m_transform;

    AABB m_meshAABB;
    mg::Real m_meshBoundingRadius;
    mg::Real m_boundingRadius;

    std::vector<CollisionShape> m_collisionShapes;

};

#endif // RENDEROBJECT_H
