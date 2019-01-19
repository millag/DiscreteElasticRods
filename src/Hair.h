#pragma once

#include "ElasticRod.h"
#include "RenderObject.h"
#include "VoxelGrid.h"
#include "HairState.h"

struct HairParams
{
	HairParams();
	~HairParams();

	mg::Real m_length;
	mg::Real m_lengthVariance;
	mg::Real m_helicalRadius;
	mg::Real m_helicalPitch;
	mg::Real m_density;
	mg::Real m_thickness;
	unsigned m_nParticles;

	mg::Vec3D m_gravity;
	mg::Real m_drag;

	bool m_resolveCollisions;
	mg::Real m_coulombFriction;

	bool m_resolveSelfInterations;
	mg::Real m_selfInterationDist;
	mg::Real m_selfStiction;
	mg::Real m_selfRepulsion;

///     constraints are enforced using PBD(Position Based Dynamics) framework
///     the parameter controls # of PBD iterations to be performed
///     NOTE that PBD DOES NOT GUARANTEE exact enforcement but converges towards the solution
///     higher value of iterations means higher precision but is more computationally expensive
	unsigned m_pbdIter;

	ElasticRodParams* m_rodParams;
};

class Hair
{
public:

	Hair();
	~Hair();

	inline unsigned getId() const { return m_id; }
	inline void setId(unsigned id) { m_id = id; }

	void getState(HairState& o_state) const;
	void setState(const HairState& state);

	void initialize();
	void reset();
	void update(mg::Real dt);

public:

	HairParams* m_params;
	const RenderObject* m_object;
	std::vector<unsigned> m_findices;
	std::vector<unsigned> m_vindices;
	std::vector<ElasticRod*> m_strands;

private:

	void resetGrid();
	void updateRod(ElasticRod& rod, mg::Real dt) const;
	void accumulateExternalForces(const ElasticRod &rod, std::vector<mg::Vec3D>& o_forces) const;
	void accumulateExternalForcesWithSelfInterations(ElasticRod &rod, std::vector<mg::Vec3D>& o_forces) const;
	void enforceConstraints(ElasticRod& rod) const;
	void enforceConstraintsWithCollision(ElasticRod& rod) const;
	void applyCollisionConstraintsIteration(ElasticRod& rod) const;

private:

	unsigned m_id;

	AABB m_volume;
	VoxelGrid* m_grid;
};
