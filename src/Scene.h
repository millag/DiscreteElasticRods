#pragma once

#include "Hair.h"
#include "Spiral.h"

class Scene
{
	friend class SceneLoader;
public:
	using MeshList = std::vector<std::shared_ptr<Mesh>>;
	using RenderObjectList = std::vector<std::shared_ptr<RenderObject>>;
	using HairList = std::vector<std::shared_ptr<Hair>>;

	void initialize();
	void update( mg::Real dt );

	inline const AABB& getBoundingVolume() const { return m_boundingVolume; }
	inline const RenderObjectList& getRenderObjects() const { return m_renderObjects; }
	inline const MeshList& getMeshes() const { return m_meshes; }

	inline std::shared_ptr<Hair> getHairById( unsigned id ) const { return ( ( id < m_hairs.size() )? m_hairs[id] : nullptr ); }

private:
	AABB m_boundingVolume;

	RenderObjectList m_renderObjects;
	MeshList m_meshes;
	HairList m_hairs;

	std::unique_ptr<Spiral> m_spiral;
};
