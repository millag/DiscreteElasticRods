#include "HairGenerator.h"
#include "Utils.h"
#include <QElapsedTimer>

void HairGenerator::generateHair( const RenderObject& skin, const std::vector<unsigned>& findices, Hair& o_hair )
{
	const auto mesh = skin.getMesh();
	assert( mesh != nullptr );

	o_hair.reset();
	o_hair.m_findices = findices;
	o_hair.m_vindices.reserve( findices.size() );

	std::set<unsigned> vIdxSet;
	for ( auto i = 0u; i < findices.size(); ++i)
	{
		const unsigned idx = mesh->getPrimitiveOffset( findices[i] );
		for ( auto j = 0u; j < mesh->getNVerticesPerPrimitive(); ++j )
		{
			const auto hasIdx = vIdxSet.insert( mesh->m_vindices[idx + j] );
			if ( hasIdx.second )
			{
				o_hair.m_vindices.push_back( *hasIdx.first );
			}
		}
	}

	o_hair.m_strands.resize( o_hair.m_vindices.size() );

	const auto normalTransform = mg::inverse( skin.getTransform() ).transpose();

	OMP_PARALLEL_LOOP
	for ( auto i = 0ll; i < static_cast<long long>( o_hair.m_vindices.size() ); ++i )
	{
		const auto idx = o_hair.m_vindices[i];
		const auto p = mg::transform_point( skin.getTransform(), mesh->m_vertices[idx] );
		const auto n = mg::transform_vector( normalTransform, mesh->m_normals[ idx ] );
		auto u = mg::Oy;
		if ( std::fabs( 1.f - std::fabs( mg::dot( n, u ) ) ) < mg::ERR )
		{
			u = mg::Ox;
		}

		generateStraightRod( o_hair.m_params, p, n, u, o_hair.m_strands[i] );
	}

	o_hair.initialize( skin );
}

void HairGenerator::generateCurlyHair( const RenderObject& skin, const std::vector<unsigned>& findices, Hair& o_hair )
{
	const auto mesh = skin.getMesh();
	assert( mesh != nullptr );

	o_hair.reset();
	o_hair.m_findices = findices;
	o_hair.m_vindices.reserve( findices.size() );

	std::set<unsigned> vIdxSet;
	for ( auto i = 0u; i < findices.size(); ++i )
	{
		const auto idx = mesh->getPrimitiveOffset(findices[i]);
		for ( auto j = 0u; j < mesh->getNVerticesPerPrimitive(); ++j )
		{
			const auto hasIdx = vIdxSet.insert( mesh->m_vindices[idx + j] );
			if ( hasIdx.second )
			{
				o_hair.m_vindices.push_back( *hasIdx.first );
			}
		}
	}

	const mg::Matrix4D normalTransform = mg::inverse( skin.getTransform() ).transpose();

	o_hair.m_strands.resize( o_hair.m_vindices.size() );

	OMP_PARALLEL_LOOP
	for ( auto i = 0ll; i < static_cast<long long>( o_hair.m_vindices.size() ); ++i )
	{
		const auto idx = o_hair.m_vindices[i];
		const auto p = mg::transform_point( skin.getTransform(), mesh->m_vertices[idx] );
		const auto n = mg::transform_vector( normalTransform, mesh->m_normals[idx] );
		auto u = mg::Oy;
		if ( std::fabs( 1.f - std::fabs( mg::dot( n, u ) ) ) < mg::ERR )
		{
			u = mg::Ox;
		}

		generateHelicalRod( o_hair.m_params, p, n, u, o_hair.m_strands[i] );
	}

	o_hair.initialize( skin );
}

void HairGenerator::generateHelicalRod( const HairParams& params,
                                        const mg::Vec3D& root,
                                        const mg::Vec3D& n,
                                        const mg::Vec3D& up,
                                        ElasticRod& o_rod )
{
	o_rod.m_ppos.resize( params.m_nParticles );
	o_rod.m_pvel.resize( params.m_nParticles );
	o_rod.m_pmass.resize( params.m_nParticles );
	o_rod.m_theta.set_size( params.m_nParticles - 1 );
	o_rod.m_isClamped.clear();

	const auto dirv = mg::normalize( mg::cross( n, up ) );
	const auto dirn = mg::normalize( n );
	const auto diru = mg::normalize( mg::cross( dirv, dirn ) );

	const auto length = params.m_length + mg::randf( -params.m_lengthVariance, params.m_lengthVariance );
	const auto volume = length * params.m_thickness * params.m_thickness * mg::Constants::pi();
	const auto pmass = params.m_density * volume / (params.m_nParticles - 1);

	const auto sign = ( mg::randf() < 0.5f )? -1.f : 1.f;
	mg::Real angle = params.m_length / std::sqrt( params.m_helicalRadius * params.m_helicalRadius + params.m_helicalPitch * params.m_helicalPitch );
	angle /= (params.m_nParticles - 1);
	const auto  phase = mg::randf( 0.f, mg::Constants::two_pi() );
	for ( auto i = 0u; i < o_rod.m_ppos.size(); ++i )
	{
//		calc point on unit circle
		o_rod.m_ppos[i] = std::cos( phase + sign * i * angle ) * diru + std::sin( phase + i * sign * angle ) * dirv
		                - std::cos( phase ) *                    diru - std::sin( phase )                    * dirv;
//		place point on helix
		o_rod.m_ppos[i] = o_rod.m_ppos[i]* params.m_helicalRadius + params.m_helicalPitch * ( i * angle ) * dirn;
//		move helix to position
		o_rod.m_ppos[i] += root;
//		init velocity and mass
		o_rod.m_pvel[i].zero();
		o_rod.m_pmass[i] = pmass;
		if ( i < static_cast<unsigned>( o_rod.m_theta.size() ) )
		{
			o_rod.m_theta( i ) = 0;
		}
	}
	o_rod.m_isClamped.insert( 0 );

	const auto e0 = o_rod.m_ppos[1] - o_rod.m_ppos[0];
	o_rod.m_u0 = mg::normalize( mg::cross( mg::cross( e0, up ), e0 ) );

	o_rod.initialize( params.m_rodParams,
	                  o_rod.m_ppos,
	                  o_rod.m_u0,
	                  o_rod.m_ppos,
	                  o_rod.m_pvel,
	                  o_rod.m_pmass,
	                  o_rod.m_theta,
	                  o_rod.m_isClamped );
}

void HairGenerator::generateStraightRod( const HairParams& params,
                                         const mg::Vec3D& root,
                                         const mg::Vec3D& n,
                                         const mg::Vec3D& up,
                                         ElasticRod& o_rod )
{
	o_rod.m_ppos.resize( params.m_nParticles );
	o_rod.m_pvel.resize( params.m_nParticles );
	o_rod.m_pmass.resize( params.m_nParticles );
	o_rod.m_theta.set_size( params.m_nParticles - 1 );
	o_rod.m_isClamped.clear();

	const auto length = params.m_length + mg::random_real( -params.m_lengthVariance, params.m_lengthVariance );
	const auto volume = length * params.m_thickness * params.m_thickness * mg::Constants::pi();
	const auto pmass = params.m_density * volume / (params.m_nParticles - 1);

//	generate straight line
	mg::Vec3D end = root + length * mg::normalize( n );
	mg::Real t = 0.0;
	for (unsigned i = 0; i < o_rod.m_ppos.size(); ++i)
	{
		t = (mg::Real)(i) / (params.m_nParticles - 1);
		o_rod.m_ppos[i] = (1 - t) * root + t * end;
//		velocity
		o_rod.m_pvel[i].zero();
//		mass
		o_rod.m_pmass[i] = pmass;
//		twist
		if (i < (unsigned)o_rod.m_theta.size())
		{
			o_rod.m_theta(i) = 0;
		}
	}
	o_rod.m_isClamped.insert(0);

	mg::Vec3D e0 = o_rod.m_ppos[1] - o_rod.m_ppos[0];
	o_rod.m_u0 = mg::cross(e0, up);
	o_rod.m_u0 = mg::normalize( mg::cross( o_rod.m_u0, e0 ) );

	o_rod.initialize( params.m_rodParams,
	                  o_rod.m_ppos,
	                  o_rod.m_u0,
	                  o_rod.m_ppos,
	                  o_rod.m_pvel,
	                  o_rod.m_pmass,
	                  o_rod.m_theta,
	                  o_rod.m_isClamped );
}
