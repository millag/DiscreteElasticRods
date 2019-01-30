#include "Hair.h"
#include "Utils.h"

Hair::Hair() = default;
Hair::~Hair()
{
	reset();
}

void Hair::getState(HairState& o_state) const
{
	o_state.m_strands.resize(m_strands.size());
	for (unsigned i = 0; i < m_strands.size(); ++i)
	{
		m_strands[i]->getState(o_state.m_strands[i]);
	}
}

void Hair::setState(const HairState& state)
{
	assert( state.m_strands.size() == m_strands.size() );
	for (unsigned i = 0; i < m_strands.size(); ++i)
	{
		m_strands[i]->setState(state.m_strands[i]);
	}
}

void Hair::reset()
{
	for ( auto it = m_strands.begin(); it != m_strands.end(); ++it )
	{
		delete (*it);
	}
	m_strands.clear();
	m_findices.clear();
	m_vindices.clear();
	m_object = nullptr;
	m_densityGrid.reset( AABB(), static_cast<VoxelGridR::size_type>( 1 ) );
	m_velocityGrid.reset( AABB(), static_cast<VoxelGridVec3D::size_type>( 1 ) );
}

void Hair::initialize()
{
	assert( m_object != nullptr );

	const auto gridCenter = m_object->getCenter();
	const auto offset = mg::Vec3D(1, 1, 1) * (m_object->getBoundingRadius() + m_params.getMaxLength());
	const AABB volume( gridCenter - offset, gridCenter + offset );
	m_densityGrid.reset( volume, m_params.m_selfInterationDist );
	m_densityGrid.clear( 0 );
	m_velocityGrid.reset( volume, m_params.m_selfInterationDist );
	m_velocityGrid.clear( mg::Vec3D(0,0,0) );
}

void Hair::updateGrid()
{
	m_densityGrid.clear( 0 );
	for ( auto it = m_strands.begin(); it != m_strands.end(); ++it )
	{
		const auto rod = *it;
		for (auto i = 1u; i < rod->m_ppos.size(); ++i)
		{
			m_densityGrid.insertValue( rod->m_ppos[i], 1 );
		}
	}

	m_velocityGrid.clear( mg::Vec3D(0,0,0) );
	for ( auto it = m_strands.begin(); it != m_strands.end(); ++it )
	{
		const auto rod = *it;
		for ( auto i = 1u; i < rod->m_ppos.size(); ++i )
		{
			m_velocityGrid.insertValue( rod->m_ppos[i], rod->m_pvel[i] );
		}
	}
}

void Hair::update(mg::Real dt)
{
	assert( m_object != nullptr );

	if (m_params.m_resolveSelfInterations)
	{
		updateGrid();
	}

	auto mesh = m_object->getMesh();

	OMP_PARALLEL_LOOP
	for ( auto i = 0ll; i < static_cast<long long>( m_vindices.size() ); ++i )
	{
		m_strands[i]->m_ppos[0] = mg::transform_point( m_object->getTransform(), mesh->m_vertices[m_vindices[i]] );
		updateRod( *m_strands[i], dt );
	}
}

// semi-implicit Euler with Verlet scheme for velocity update
void Hair::updateRod(ElasticRod& rod, mg::Real dt) const
{
	std::vector<mg::Vec3D> forces( rod.m_ppos.size(), mg::Vec3D(0,0,0) );
	rod.accumulateInternalElasticForces(forces);

	if (m_params.m_resolveSelfInterations)
	{
		accumulateExternalForcesWithSelfInterations(rod, forces);
	} else
	{
		accumulateExternalForces(rod, forces);
	}

//    integrate centerline - semi- implicit Euler
	std::vector<mg::Vec3D> prevPos(rod.m_ppos.size());
	for (unsigned i = 1; i < rod.m_ppos.size(); ++i)
	{
		prevPos[i] = rod.m_ppos[i];
		rod.m_pvel[i] += dt * forces[i] / rod.m_pmass[i];
		rod.m_ppos[i] += rod.m_pvel[i] * dt;
	}

	if (m_params.m_resolveCollisions)
	{
		enforceConstraintsWithCollision(rod);
	} else
	{
		enforceConstraints(rod);
	}

//    velocity correction using Verlet scheme:
	for (unsigned i = 1; i < rod.m_ppos.size(); ++i)
	{
		rod.m_pvel[i] = (rod.m_ppos[i] - prevPos[i]) / dt;
	}

	rod.updateCurrentState();
}

void Hair::enforceConstraints(ElasticRod& rod) const
{
	for (auto k = 0u; k < m_params.m_pbdIter; ++k)
	{
		rod.applyInternalConstraintsIteration();
	}
}

void Hair::enforceConstraintsWithCollision(ElasticRod& rod) const
{
	for (auto k = 0u; k < m_params.m_pbdIter; ++k)
	{
		applyCollisionConstraintsIteration(rod);
		rod.applyInternalConstraintsIteration();
	}
}

void Hair::applyCollisionConstraintsIteration(ElasticRod& rod) const
{
	mg::Vec3D collision_p, normal;
	for ( unsigned i = 1; i < rod.m_ppos.size(); ++i)
	{
		if (m_object->isInsideObject(rod.m_ppos[i], collision_p, normal))
		{
			rod.m_ppos[i] = collision_p;
		}
	}
}

void Hair::accumulateExternalForcesWithSelfInterations(ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
	const auto dr = m_densityGrid.getVolexSize() * 0.5;
	const auto dr_x2 = m_densityGrid.getVolexSize();

	mg::Vec3D p1, p2, forceDensity;
	for ( auto i = 1u; i < o_forces.size(); ++i )
	{
//		gravity
		o_forces[i] += m_params.m_gravity * rod.m_pmass[i];
//		drag
		o_forces[i] -= m_params.m_drag * rod.m_pvel[i].length() * rod.m_pvel[i];

//		self interactions:
		if( m_params.m_selfStiction > mg::ERR )
		{
//			self stiction acts as averaging of the velocity with nearby particles
//			the velocity however is not direcly averaged by trilinear interpolation
//			is used to calculate the avg. value for particle's current position
			const auto density = m_densityGrid.estimateValueAt( rod.m_ppos[i] ) ;
			if ( density > mg::ERR )
			{
				const mg::Vec3D v = m_velocityGrid.estimateValueAt( rod.m_ppos[i] ) / density;
				rod.m_pvel[i] = mg::lerp( rod.m_pvel[i], v, m_params.m_selfStiction );
			}
		}

		if( m_params.m_selfRepulsion > mg::ERR )
		{
//			self repulsion
//			NOTE: forceDensity = -pressureGradient, where pressureGradient is estimated
//			numerically from the density grid
			p1.set( rod.m_ppos[i][0] - dr, rod.m_ppos[i][1], rod.m_ppos[i][2] );
			p2.set( rod.m_ppos[i][0] + dr, rod.m_ppos[i][1], rod.m_ppos[i][2] );
			forceDensity[0] = (m_densityGrid.estimateValueAt( p1 ) - m_densityGrid.estimateValueAt( p2 )) / dr_x2;

			p1.set( rod.m_ppos[i][0], rod.m_ppos[i][1] - dr, rod.m_ppos[i][2] );
			p2.set( rod.m_ppos[i][0], rod.m_ppos[i][1] + dr, rod.m_ppos[i][2] );
			forceDensity[1] = (m_densityGrid.estimateValueAt( p1 ) - m_densityGrid.estimateValueAt( p2 )) / dr_x2;

			p1.set( rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] - dr );
			p2.set( rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] + dr );
			forceDensity[2] = (m_densityGrid.estimateValueAt( p1 ) - m_densityGrid.estimateValueAt( p2 )) / dr_x2;

			o_forces[i] += m_params.m_selfRepulsion * forceDensity;
		}
	}
}

void Hair::accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
	for (unsigned i = 1; i < o_forces.size(); ++i)
	{
//        gravity
		o_forces[i] += m_params.m_gravity * rod.m_pmass[i];
//        drag
		o_forces[i] -= m_params.m_drag * rod.m_pvel[i].length() * rod.m_pvel[i];
	}
}




//for (SIter it = m_strands.begin(); it != m_strands.end(); ++it)
//{
//    Strand& strand = *it;
//    for (unsigned i = 1; i < strand.m_ppos.size(); ++i)
//    {
////            calculate self friction
//        strand.m_pvel[i] =  (1 - m_selfFriction) * strand.m_pvel[i] +
//                            m_selfFriction * m_grid->getInterpolatedVelocity(strand.m_ppos[i]) / m_grid->getInterpolatedDensity(strand.m_ppos[i]);

////            calculate self repulsion
//        p1.set(strand.m_ppos[i][0] - dr, strand.m_ppos[i][1], strand.m_ppos[i][2]);
//        p2.set(strand.m_ppos[i][0] + dr, strand.m_ppos[i][1], strand.m_ppos[i][2]);
//        pressureGradient[0] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        p1.set(strand.m_ppos[i][0], strand.m_ppos[i][1] - dr, strand.m_ppos[i][2]);
//        p2.set(strand.m_ppos[i][0], strand.m_ppos[i][1] + dr, strand.m_ppos[i][2]);
//        pressureGradient[1] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        p1.set(strand.m_ppos[i][0], strand.m_ppos[i][1], strand.m_ppos[i][2] - dr);
//        p2.set(strand.m_ppos[i][0], strand.m_ppos[i][1], strand.m_ppos[i][2] + dr);
//        pressureGradient[2] = (m_grid->getInterpolatedDensity(p1) - m_grid->getInterpolatedDensity(p2)) * 0.5 / dr;

//        strand.m_pvel[i] = 0.5 * (2 * strand.m_pvel[i] +  m_selfRepulsion * dt * pressureGradient / strand.m_pmass[i]);
//    }
//}

//        add friction
//        if (rod.m_isClamped.count(i))
//        {
//            continue;
//        }
//        mg::Vec3D collision_p, normal, vt, vn;
//        mg::Real dist = m_object->distanceToSurface(rod.m_ppos[i], collision_p, normal);
//        if (dist > m_params->m_thickness + 0.0001)
//        {
//            continue;
//        }

//        vn = mg::dot(normal, rod.m_pvel[i]) * normal;
//        vt = rod.m_pvel[i] + vn;

//        if (vt.length_squared() > mg::ERR)
//        {
//            rod.m_pvel[i] -= std::max(1 - m_params->m_coulombFriction * vn.length() / vt.length(), (mg::Real)0) * vt;
//        }
