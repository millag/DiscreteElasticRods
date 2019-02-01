#include "Spiral.h"
#include "Utils.h"

void Spiral::initialize( const RenderObject& object )
{
	m_object = &object;

	m_rodParams1 = ElasticRodParams( 0.4, 0.7, 1000, ElasticRodParams::MinimizationStrategy::Newton );
	m_rodParams2 = ElasticRodParams( 0.4, 0.7, 1000, ElasticRodParams::MinimizationStrategy::BFGS );
	m_rodParams3 = ElasticRodParams( 0.4, 0.7, 1000, ElasticRodParams::MinimizationStrategy::None );

	std::vector<mg::Vec3D> restpos( m_params.m_nParticles );
	std::vector<mg::Vec3D> pos( m_params.m_nParticles );
	std::vector<mg::Vec3D> vel( m_params.m_nParticles );
	std::vector<mg::Real> mass( m_params.m_nParticles );
	ColumnVector theta( m_params.m_nParticles - 1);
	std::set<unsigned> isClamped;

	const auto center =  m_object->getCenter();
	const auto dirx = mg::Ox;
	const auto diry = -mg::Oy;
	const auto dirz = mg::Oz;

	const auto angle = m_params.m_lenght / ((m_params.m_nParticles - 1) * m_params.m_radius);
	const auto step = 0.1f;
	auto len = 0.f;
	for (auto i = 0u; i < restpos.size(); ++i)
	{
		restpos[i] = (std::cos( i * angle ) * dirx  + std::sin( i * angle ) * dirz - dirx) * m_params.m_radius + center;
		restpos[i] += diry * i * step;
		if ( i > 0 )
		{
			len += mg::length( restpos[i] - restpos[i - 1] );
		}
	}
	isClamped.insert( 0 );

	const auto e0 = restpos[1] - restpos[0];
	const auto u0 = mg::normalize( mg::cross( mg::cross( e0, dirx ), e0 ) );

	const auto start = restpos[0];
	const mg::Vec3D end = start + len * dirx;
	for ( auto i = 0u; i < pos.size(); ++i )
	{
		const auto t = static_cast<mg::Real>( i ) / (m_params.m_nParticles - 1);
		pos[i] = mg::lerp( start, end, t );

		vel[i].zero();
		mass[i] = 0.05;
		if ( i < static_cast<unsigned>( theta.size() ) )
		{
			theta( i ) = 0;
		}
	}

	m_strands.resize( 3 );

	m_strands[0].initialize( m_rodParams1, restpos, u0, pos, vel, mass, theta, isClamped );

	for ( auto i = 0u; i < pos.size(); ++i )
	{
		pos[i] -= m_params.m_offset * mg::Ox;
	}
	m_strands[1].initialize( m_rodParams1, restpos, u0, pos, vel, mass, theta, isClamped );

	for ( auto i = 0u; i < pos.size(); ++i )
	{
		pos[i] += 2 * m_params.m_offset * mg::Ox;
	}
	m_strands[2].initialize(m_rodParams3, restpos, u0, pos, vel, mass, theta, isClamped );
}

void Spiral::update(mg::Real dt)
{
	const auto center = m_object->getCenter();

	m_strands[0].m_ppos[0] = center;
	updateRod( m_strands[0], dt );

	m_strands[1].m_ppos[0] = center - m_params.m_offset * mg::Ox;
	updateRod( m_strands[1], dt );

	m_strands[2].m_ppos[0] = center + m_params.m_offset * mg::Ox;
	updateRod( m_strands[2], dt );
}

/* semi-implicit Euler with Verlet scheme for velocity update */
void Spiral::updateRod( ElasticRod& rod, mg::Real dt ) const
{
	std::vector<mg::Vec3D> forces( rod.m_ppos.size(), mg::Vec3D(0,0,0) );
	rod.accumulateInternalElasticForces( forces );
	accumulateExternalForces( rod, forces );

//	integrate centerline - semi- implicit Euler
	std::vector<mg::Vec3D> prevPos( rod.m_ppos.size() );
	for ( auto i = 0u; i < rod.m_ppos.size(); ++i )
	{
		prevPos[i] = rod.m_ppos[i];
		rod.m_pvel[i] += dt * forces[i] / rod.m_pmass[i];
		rod.m_ppos[i] += dt * rod.m_pvel[i];
	}

	for ( auto k = 0u; k < m_params.m_pbdIter; ++k )
	{
		rod.applyInternalConstraintsIteration();
	}

//	velocity correction using Verlet scheme
	for ( auto i = 0u; i < rod.m_ppos.size(); ++i )
	{
		rod.m_pvel[i] = (rod.m_ppos[i] - prevPos[i]) / dt;
	}

	rod.updateCurrentState();
}

void Spiral::accumulateExternalForces( const ElasticRod& rod, std::vector<mg::Vec3D>& o_forces ) const
{
	for ( auto i = 0u; i < o_forces.size(); ++i )
	{
		if ( !rod.m_isClamped.count( i ) )
		{
			o_forces[i] += mg::Gravity * rod.m_pmass[i] - 0.001 * rod.m_pvel[i];
		}
	}
}

//void Spiral::initialize(const RenderObject& object)
//{
//    m_object = &object;

//    std::vector<mg::Vec3D> restpos(m_nParticles);
//    std::vector<mg::Vec3D> pos(m_nParticles);
//    std::vector<mg::Vec3D> vel(m_nParticles);
//    std::vector<mg::Real> mass(m_nParticles);
//    ColumnVector theta(m_nParticles - 1);
//    std::set<unsigned> isClamped;

//    mg::Vec3D center(m_object->getCenter());
//    mg::Vec3D dirx = mg::Vec3D(1, 0, 0);
//    mg::Vec3D diry = mg::Vec3D(0, 1, 0);
//    mg::Vec3D dirz = mg::Vec3D(0, 0, 1);

//    mg::Real radius = m_lenght / mg::Constants::pi();
//    mg::Real angle = mg::Constants::pi() / (m_nParticles - 1);
//    for (unsigned i = 0; i < restpos.size(); ++i)
//    {
//        unsigned j = i;//(m_nParticles - 1) - i;
//        pos[i] = (std::cos(j * angle) * dirx  + std::sin(j * angle) * diry - dirx) * radius + center;
//    }
//    isClamped.insert(0);
////    isClamped.insert(m_nParticles - 1);

//    mg::Vec3D start = pos[0];
//    mg::Vec3D end = start + m_lenght * dirx;
//    mg::Real t = 0.0;
//    for (unsigned i = 0; i < pos.size(); ++i)
//    {
//        t = (mg::Real)(i) / (m_nParticles - 1);
//        restpos[i] = (1 - t) * start + t * end;

//        vel[i].zero();
//        mass[i] = 0.001;
//        if (i < (unsigned)theta.size())
//        {
//            theta(i) = 0;
//        }
//    }

//    m_strands.resize(1);
//    for (unsigned i = 0; i < m_strands.size(); ++i)
//    {
////        m_strands[i].init(pos, diry, pos, vel, mass, theta, isClamped);
//        m_strands[i].init(restpos, -dirz, pos, vel, mass, theta, isClamped);
//    }
//}
