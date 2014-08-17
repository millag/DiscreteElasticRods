#ifndef SCENE_H
#define SCENE_H

#include <string>
#include <vector>
#include <string>
#include <map>

#include "IServant.h"
#include "HairStrand.h"
#include "Spiral.h"

class Scene : virtual public IServant {
public:
    Scene();
    ~Scene();

    void initialize();
    void update(mg::Real dt);

    const AABB& getBoundingVolume() const { return m_boundingVolume; }
    const std::vector<ElasticRod*>& getStrands() const { return m_spiral->getStrands(); }
    const std::vector<RenderObject*>& getRenderObjects() const { return m_renderObjects; }
    const std::vector<Mesh*>& getMeshes() const { return m_meshes; }
    void findObjectsWithinDistance(const mg::Vec3D& pos, mg::Real dist, std::vector<RenderObject*>& o_objects);

//    const Grid& getGrid() const;

private:
    AABB m_boundingVolume;

    HairStrand* m_hair;
    Spiral* m_spiral;
    std::vector<RenderObject*> m_renderObjects;
    std::vector<Mesh*> m_meshes;
    std::map<std::string, unsigned> m_meshMap;
};

#endif // SCENE_H
