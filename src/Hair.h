#pragma once

#include "ElasticRod.h"
#include "RenderObject.h"
#include "VoxelGrid.h"

struct HairParams
{
	mg::Real getMaxLength() const
	{
		return m_length + m_lengthVariance;
	}

	mg::Real m_length = 6.f;
	mg::Real m_lengthVariance = 1.f;
	mg::Real m_helicalRadius = 0.3f;
	mg::Real m_helicalPitch = 0.15f;
	mg::Real m_density = 0.001f;
	mg::Real m_thickness = 0.07f;
	unsigned m_nParticles = 15;

	mg::Vec3D m_gravity = mg::Gravity;
	mg::Real m_drag = 0.0000003f;

	bool m_resolveCollisions = true;
	mg::Real m_coulombFriction = 0.2f;

	bool m_resolveSelfInterations = true;
	mg::Real m_selfInterationDist = 0.4f;
	mg::Real m_selfStiction = 0.2f;
	mg::Real m_selfRepulsion = 0.000005f;

///     constraints are enforced using PBD(Position Based Dynamics) framework
///     the parameter controls # of PBD iterations to be performed
///     NOTE that PBD DOES NOT GUARANTEE exact enforcement but converges towards the solution
///     higher value of iterations means higher precision but is more computationally expensive
	unsigned m_pbdIter = 6;

	ElasticRodParams m_rodParams = ElasticRodParams( 0.00006, 0.0005, 1000, ElasticRodParams::MinimizationStrategy::NONE );
};

struct HairState
{
	friend class Hair;
public:
	inline void clear()
	{
		m_strands.clear();
	}

private:
	std::vector<ElasticRodState> m_strands;
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
	HairParams m_params;
	const RenderObject* m_object = nullptr;
	std::vector<unsigned> m_findices;
	std::vector<unsigned> m_vindices;
	std::vector<ElasticRod*> m_strands;

private:
	unsigned m_id = 0;
	VoxelGrid m_grid;

private:
	void updateGrid();
	void updateRod(ElasticRod& rod, mg::Real dt) const;
	void accumulateExternalForces(const ElasticRod &rod, std::vector<mg::Vec3D>& o_forces) const;
	void accumulateExternalForcesWithSelfInterations(ElasticRod &rod, std::vector<mg::Vec3D>& o_forces) const;
	void enforceConstraints(ElasticRod& rod) const;
	void enforceConstraintsWithCollision(ElasticRod& rod) const;
	void applyCollisionConstraintsIteration(ElasticRod& rod) const;
};
