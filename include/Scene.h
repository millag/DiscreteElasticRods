#ifndef SCENE_H
#define SCENE_H

#include <string>
#include <vector>
#include <string>
#include <map>

#include "Hair.h"
#include "Spiral.h"

class Loader;

class Scene
{
public:
    Scene();
    ~Scene();

    void initialize();
    void update(mg::Real dt);

    inline const AABB& getBoundingVolume() const { return m_boundingVolume; }
    inline const std::vector<ElasticRod*>& getStrands() const { return m_hair->m_strands; }
    inline const std::vector<RenderObject*>& getRenderObjects() const { return m_renderObjects; }
    inline const std::vector<Mesh*>& getMeshes() const { return m_meshes; }

    void findObjectsWithinDistance(const mg::Vec3D& pos, mg::Real dist, std::vector<RenderObject*>& o_objects);

private:
    AABB m_boundingVolume;

    Hair* m_hair;
    Spiral* m_spiral;
    std::vector<RenderObject*> m_renderObjects;
    std::vector<Mesh*> m_meshes;

private:
    friend class Loader;
};

#endif // SCENE_H
