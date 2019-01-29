#include "VoxelGrid.h"
#include <array>

void VoxelGrid::reset(const AABB &volume, size_type resolution )
{
	assert( resolution > 0 );

	m_volume = volume;
	m_resolution = resolution + 1;
	m_voxelSize = m_volume.getWidth() / ( m_resolution - 1 );
	m_voxels.resize( m_resolution * m_resolution * m_resolution );

	clear();
}

void VoxelGrid::reset(const AABB &volume, mg::Real voxelSize )
{
	const auto w = volume.getWidth();
	voxelSize = std::min( std::fabs( voxelSize ), w );

	m_volume = volume;
	m_resolution = static_cast<size_type>( w / voxelSize ) + 1;
	m_voxelSize = w / ( m_resolution - 1 );
	m_voxels.resize( m_resolution * m_resolution * m_resolution );

	clear();
}

void VoxelGrid::clear()
{
	OMP_PARALLEL_LOOP
	for ( auto i = 0ll; i < static_cast<long long>( m_voxels.size() ); ++i )
	{
		m_voxels[i].m_density = 0;
		m_voxels[i].m_velocity.zero();
	}
}

void VoxelGrid::insertDensity( const mg::Vec3D& pos, mg::Real val )
{
	size_type i,j,k;
	mg::Vec3D localPos;
	findVoxel( pos, i, j, k, localPos );

	assert( getVoxelIdx( i, j, k ) < m_voxels.size() );

	m_voxels[getVoxelIdx( i  , j  , k   )].m_density += (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j  , k   )].m_density += (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j+1, k   )].m_density += (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j+1, k   )].m_density += (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j  , k+1 )].m_density +=      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j  , k+1 )].m_density +=      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j+1, k+1 )].m_density +=      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j+1, k+1 )].m_density +=      localPos[2]  *      localPos[1]  *      localPos[0]  * val;
}

void VoxelGrid::insertVelocity( const mg::Vec3D& pos, const mg::Vec3D& val )
{
	size_type i, j, k;
	mg::Vec3D localPos;
	findVoxel(pos, i, j, k, localPos);

	assert( getVoxelIdx( i, j, k ) < m_voxels.size() );

	m_voxels[getVoxelIdx( i  , j  , k   )].m_velocity += (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j  , k   )].m_velocity += (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j+1, k   )].m_velocity += (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j+1, k   )].m_velocity += (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j  , k+1 )].m_velocity +=      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j  , k+1 )].m_velocity +=      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * val;
	m_voxels[getVoxelIdx( i  , j+1, k+1 )].m_velocity +=      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * val;
	m_voxels[getVoxelIdx( i+1, j+1, k+1 )].m_velocity +=      localPos[2]  *      localPos[1]  *      localPos[0]  * val;
}

void VoxelGrid::getInterpolatedDensity(const mg::Vec3D& pos, mg::Real& o_density) const
{
	size_type i,j,k;
	mg::Vec3D localPos;
	findVoxel(pos, i, j, k, localPos);

	assert( getVoxelIdx( i, j, k ) < m_voxels.size() );

	o_density = 0;
	o_density += (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k )].m_density;
	o_density += (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k )].m_density;
	o_density += (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k )].m_density;
	o_density += (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k )].m_density;
	o_density +=      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k+1 )].m_density;
	o_density +=      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k+1 )].m_density;
	o_density +=      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k+1 )].m_density;
	o_density +=      localPos[2]  *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k+1 )].m_density;
}

void VoxelGrid::getInterpolatedVelocity(const mg::Vec3D& pos, mg::Vec3D& o_velocity) const
{
	size_type i,j,k;
	mg::Vec3D localPos;
	findVoxel(pos, i, j, k, localPos);

	assert( getVoxelIdx( i, j, k ) < m_voxels.size() );

	o_velocity.zero();
	o_velocity += (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k )].m_velocity;
	o_velocity += (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k )].m_velocity;
	o_velocity += (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k )].m_velocity;
	o_velocity += (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k )].m_velocity;
	o_velocity +=      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k+1 )].m_velocity;
	o_velocity +=      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k+1 )].m_velocity;
	o_velocity +=      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k+1 )].m_velocity;
	o_velocity +=      localPos[2]  *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k+1 )].m_velocity;
}

void VoxelGrid::findVoxel( const mg::Vec3D& pos, size_type& o_i, size_type& o_j, size_type& o_k, mg::Vec3D& o_localPos ) const
{
	o_localPos = pos - m_volume.getMin();

	constexpr const auto Dim = 3;
	std::array<mg::Real, Dim> voxelIdx;
	for ( auto i = 0; i < Dim; ++i )
	{
		voxelIdx[i] = std::floor( o_localPos[i] / m_voxelSize );
	}

	const auto maxIdx = static_cast<mg::Real>( m_resolution - 1 );
	for ( auto i = 0; i < Dim; ++i )
	{
		o_localPos[i] = mg::clamp( std::fmod( o_localPos[i], m_voxelSize ) / m_voxelSize, 0.f, 1.f );

		if ( voxelIdx[i] < 0.f )
		{
			voxelIdx[i] = 0.f;
			o_localPos[i] = 0.f;
		}
		else if ( voxelIdx[i] >= maxIdx )
		{
			voxelIdx[i] = maxIdx - 1.f;
			o_localPos[i] = 1.f;
		}
	}

	o_i = static_cast< size_type >( voxelIdx[0] );
	o_j = static_cast< size_type >( voxelIdx[1] );
	o_k = static_cast< size_type >( voxelIdx[2] );
}
