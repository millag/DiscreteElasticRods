#pragma once

#include "Hair.h"
#include "Spiral.h"

class Scene
{
	friend class SceneLoader;
public:
	Scene() = default;
	~Scene();

	void initialize();
	void update(mg::Real dt);

	inline const AABB& getBoundingVolume() const { return m_boundingVolume; }
	inline const std::vector<RenderObject*>& getRenderObjects() const { return m_renderObjects; }
	inline const std::vector<Mesh*>& getMeshes() const { return m_meshes; }

	inline Hair* getHairById(unsigned id) const { return ( ( id < m_hairs.size() )? m_hairs[id] : nullptr ); }

private:
	AABB m_boundingVolume;

	std::vector<Hair*> m_hairs;
	std::vector<RenderObject*> m_renderObjects;
	std::vector<Mesh*> m_meshes;

	std::unique_ptr<Spiral> m_spiral;
};
