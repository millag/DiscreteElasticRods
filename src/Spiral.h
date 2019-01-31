#pragma once

#include "RenderObject.h"
#include "ElasticRod.h"

struct SpiralParams
{
	mg::Real m_radius = 0.2f;
	mg::Real m_lenght = 4.f;
	mg::Real m_offset = 3.f;
	unsigned m_nParticles = 20;
	unsigned m_pbdIter = 4;
};

class Spiral
{
public:
	void initialize( const RenderObject& object );
	void update( mg::Real dt );

public:
	SpiralParams m_params;

private:
	const RenderObject* m_object = nullptr;
	ElasticRodParams m_rodParams1;
	ElasticRodParams m_rodParams2;
	ElasticRodParams m_rodParams3;
	std::vector<ElasticRod> m_strands;

private:
	void updateRod(ElasticRod& rod, mg::Real dt) const;
	void accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const;
};
