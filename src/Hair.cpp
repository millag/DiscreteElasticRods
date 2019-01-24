#include "Hair.h"
#include "Utils.h"

Hair::Hair(): m_object(nullptr), m_id(-1), m_grid(nullptr)
{
//    init default hair params
	m_params = new HairParams();
	m_params->m_length = 6;
	m_params->m_lengthVariance = 1;
	m_params->m_helicalRadius = 0.3;
	m_params->m_helicalPitch = 0.15;
	m_params->m_density = 0.001;
	m_params->m_thickness = 0.07;
	m_params->m_nParticles = 15;

	m_params->m_gravity.set(0, -9.81, 0);
	m_params->m_drag = 0.0000003;

	m_params->m_resolveCollisions = 1;
	m_params->m_coulombFriction = 0.2;

	m_params->m_resolveSelfInterations = 1;
	m_params->m_selfInterationDist = 0.4;
	m_params->m_selfStiction = 0.2;
	m_params->m_selfRepulsion = 0.000005;

	m_params->m_pbdIter = 6;

	mg::Real bendStiffness = 0.00006;
	mg::Real twistStiffness = 0.0005;
	mg::Real maxElasticForce = 1000;
	m_params->m_rodParams = ElasticRodParams( bendStiffness, twistStiffness, maxElasticForce, ElasticRodParams::NONE );
}

Hair::~Hair()
{
	reset();
	delete m_params;
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
	if (m_grid != nullptr)
	{
		delete m_grid;
	}
	typedef std::vector<ElasticRod*>::iterator Iter;
	for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
	{
		delete (*it);
	}
	m_strands.clear();
	m_findices.clear();
	m_vindices.clear();
	m_object = nullptr;
	m_grid = nullptr;
}

void Hair::initialize()
{
	assert( m_object != nullptr );
	assert( m_grid == nullptr );

	mg::Vec3D gridCenter = m_object->getCenter();
	mg::Vec3D offset = (m_object->getBoundingRadius() + m_params->m_length + m_params->m_lengthVariance) * mg::Vec3D(1, 1, 1);
	m_volume.reshape(gridCenter - offset, gridCenter + offset);

	m_grid = new VoxelGrid(m_volume, m_volume.getWidth() / m_params->m_selfInterationDist);
	m_grid->initialize();
}

void Hair::resetGrid()
{
	mg::Vec3D gridCenter = m_object->getCenter();
	mg::Vec3D offset = (m_object->getBoundingRadius() + m_params->m_length + m_params->m_lengthVariance) * mg::Vec3D(1, 1, 1);
	m_volume.reshape(gridCenter - offset, gridCenter + offset);

	m_grid->reset();

	typedef std::vector<ElasticRod*>::iterator Iter;
	for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
	{
		ElasticRod* rod = *it;
		for (unsigned i = 1; i < rod->m_ppos.size(); ++i)
		{
			m_grid->insertDensity(rod->m_ppos[i], 1);
		}
	}

	for (Iter it = m_strands.begin(); it != m_strands.end(); ++it)
	{
		ElasticRod* rod = *it;
		for (unsigned i = 1; i < rod->m_ppos.size(); ++i)
		{
			m_grid->insertVelocity(rod->m_ppos[i], rod->m_pvel[i]);
		}
	}

}

void Hair::update(mg::Real dt)
{
	assert( m_object != nullptr );
	assert( m_grid != nullptr );

	if (m_params->m_resolveSelfInterations)
	{
		resetGrid();
	}

	const Mesh* mesh = m_object->getMesh();

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
	std::vector<mg::Vec3D> forces(rod.m_ppos.size(), mg::Vec3D(0,0,0));
	rod.accumulateInternalElasticForces(forces);

	if (m_params->m_resolveSelfInterations)
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

	if (m_params->m_resolveCollisions)
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
	for (unsigned k = 0; k < m_params->m_pbdIter; ++k)
	{
		rod.applyInternalConstraintsIteration();
	}
}

void Hair::enforceConstraintsWithCollision(ElasticRod& rod) const
{
	for (unsigned k = 0; k < m_params->m_pbdIter; ++k)
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
	mg::Vec3D p1(0,0,0), p2(0,0,0), pressureGradient(0,0,0);
	mg::Real dr = m_params->m_selfInterationDist;
	mg::Real density_p1, density_p2;

	for (unsigned i = 1; i < o_forces.size(); ++i)
	{
//        gravity
		o_forces[i] += m_params->m_gravity * rod.m_pmass[i];
//        drag
		o_forces[i] -= m_params->m_drag * rod.m_pvel[i].length() * rod.m_pvel[i];


//        self interactions:

//        self stiction acts as averaging of the velocity for near by particles
//        the velocity however is not direcly averaged by trilinear interpolation is used to calculate the avg. value for particle's current position
		m_grid->getInterpolatedDensity(rod.m_ppos[i], density_p1);
		if ( density_p1 > mg::ERR)
		{
			m_grid->getInterpolatedVelocity(rod.m_ppos[i], p1);
			rod.m_pvel[i] =  (1 - m_params->m_selfStiction) * rod.m_pvel[i] + m_params->m_selfStiction * p1 / density_p1;
		}

//        self repulsion
		p1.set(rod.m_ppos[i][0] - dr, rod.m_ppos[i][1], rod.m_ppos[i][2]);
		p2.set(rod.m_ppos[i][0] + dr, rod.m_ppos[i][1], rod.m_ppos[i][2]);
		m_grid->getInterpolatedDensity(p1, density_p1);
		m_grid->getInterpolatedDensity(p2, density_p2);
		pressureGradient[0] = ( density_p1 - density_p2 ) * 0.5 / dr;

		p1.set(rod.m_ppos[i][0], rod.m_ppos[i][1] - dr, rod.m_ppos[i][2]);
		p2.set(rod.m_ppos[i][0], rod.m_ppos[i][1] + dr, rod.m_ppos[i][2]);
		m_grid->getInterpolatedDensity(p1, density_p1);
		m_grid->getInterpolatedDensity(p2, density_p2);
		pressureGradient[1] = ( density_p1 - density_p2 ) * 0.5 / dr;

		p1.set(rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] - dr);
		p2.set(rod.m_ppos[i][0], rod.m_ppos[i][1], rod.m_ppos[i][2] + dr);
		m_grid->getInterpolatedDensity(p1, density_p1);
		m_grid->getInterpolatedDensity(p2, density_p2);
		pressureGradient[2] = ( density_p1 - density_p2) * 0.5 / dr;

		o_forces[i] += m_params->m_selfRepulsion * pressureGradient;

	}
}

void Hair::accumulateExternalForces(const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces) const
{
	for (unsigned i = 1; i < o_forces.size(); ++i)
	{
//        gravity
		o_forces[i] += m_params->m_gravity * rod.m_pmass[i];
//        drag
		o_forces[i] -= m_params->m_drag * rod.m_pvel[i].length() * rod.m_pvel[i];
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
