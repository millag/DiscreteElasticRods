#pragma once

#include "RenderObject.h"
#include "ElasticRod.h"

class Spiral
{
public:
	Spiral();
	~Spiral();

	void initialize(const RenderObject& object);
	void update(mg::Real dt);

public:
	mg::Real m_radius;
	mg::Real m_lenght;

	mg::Real m_offset;

	unsigned m_nParticles;
	unsigned m_pbdIter;

	std::vector<ElasticRod*> m_strands;

private:
	ElasticRodParams m_rodParams1;
	ElasticRodParams m_rodParams2;
	ElasticRodParams m_rodParams3;
	const RenderObject* m_object;

private:
	void updateRod(ElasticRod& rod, mg::Real dt) const;
	void accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const;
};
