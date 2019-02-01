#include "HairGenerator.h"
#include "Utils.h"
#include <QElapsedTimer>


static void generateStraightRod( mg::Real length,
                                 const mg::Vec3D& root,
                                 const mg::Vec3D& n,
                                 const mg::Vec3D& up,
                                 ElasticRod& o_rod )
{
	UNUSED_VALUE( up );

	const auto endPos = mg::Vec3D( root + length * mg::normalize( n ) );
	for ( auto i = 0u; i < o_rod.m_ppos.size(); ++i )
	{
		const auto t = static_cast<mg::Real>( i ) / (o_rod.m_ppos.size() - 1);
		o_rod.m_ppos[i] = mg::lerp( root, endPos, t );
	}
}

static void generateHelicalRod( const HairParams& params,
                                const mg::Vec3D& root,
                                const mg::Vec3D& n,
                                const mg::Vec3D& up,
                                ElasticRod& o_rod )
{
	const auto dirv = mg::normalize( mg::cross( n, up ) );
	const auto dirn = mg::normalize( n );
	const auto diru = mg::normalize( mg::cross( dirv, dirn ) );

	const auto arcLen = mg::Constants::two_pi() * params.m_helicalRadius;
	const auto revArcLen = std::sqrt( arcLen * arcLen + params.m_helicalPitch * params.m_helicalPitch );
	const auto nTurns = params.m_length / revArcLen;
	const auto totalRevAngle = mg::Constants::two_pi() * nTurns;
	const auto sign = ( mg::randf() < 0.5f )? -1.f : 1.f;
	const auto phase = mg::randf( 0.f, mg::Constants::two_pi() );
	for ( auto i = 0u; i < o_rod.m_ppos.size(); ++i )
	{
		const auto t = static_cast<mg::Real>( i ) / (o_rod.m_ppos.size() - 1);
//		calc point on unit circle and shift it by phase
		o_rod.m_ppos[i] = std::cos( phase + sign * t * totalRevAngle ) * diru + std::sin( phase + t * sign * totalRevAngle ) * dirv
		                - std::cos( phase ) *                    diru - std::sin( phase )                    * dirv;
//		place point on helix
		o_rod.m_ppos[i] = params.m_helicalRadius * o_rod.m_ppos[i] + params.m_helicalPitch * t * nTurns * dirn;
//		move helix to position
		o_rod.m_ppos[i] += root;
	}
}

void HairGenerator::generateHair( const RenderObject& skin,
                                  const std::vector<unsigned>& findices,
                                  Hair& o_hair )
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

		generateRod( o_hair.m_params, p, n, u, o_hair.m_strands[i] );
	}

	o_hair.initialize( skin );
}

void HairGenerator::generateRod( const HairParams& params,
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

	if ( mg::almostEqual( params.m_helicalRadius, 0.f ) )
	{
		generateStraightRod( length, root, n, up, o_rod );
	}
	else
	{
		generateHelicalRod( params, root, n, up, o_rod );
	}

	for ( auto i = 0u; i < o_rod.m_pvel.size(); ++i )
	{
		o_rod.m_pvel[i].zero();
	}

	for ( auto i = 0u; i < o_rod.m_pmass.size(); ++i )
	{
		 o_rod.m_pmass[i] = pmass;
	}

	o_rod.m_theta = 0;
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
