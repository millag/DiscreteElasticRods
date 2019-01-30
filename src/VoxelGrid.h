#pragma once

#include "AABB.h"
#include <array>

/// Regularly spaced grid encompassing some volume in 3D space. Each voxel
/// represents a single lattice point on the grid and not a volume and holds
/// the sample value at that point. How many voxels/samples the grid has depends
/// on the grid resolution - f.e. a VoxelGrid with non-zero volume and resolution 1
/// has 8 samples, one per each vertex of the volume bbox. The space between each
/// voxel is not represented in the dataset, but a value at arbitrary point in 3D
/// space is approximated via tri-linear interpolation. Inserting a value in the
/// VoxelGrid it's distributed among the 8 surrounding lattice points using
/// again tri-linear interpolation
template <typename T>
class VoxelGrid
{
	using VoxelList = std::vector<T>;

public:
	using value_type = T;
	using size_type = unsigned;

	constexpr const AABB& getVolume() const { return m_volume; }
	constexpr size_type getResolution() const { return m_resolution - 1; }
	constexpr mg::Real getVolexSize() const { return m_voxelSize; }

	/// Setsup and allocate memory for grid with the given volume and resolution
	/// Note that resolution should be set to at least 1 and voxels do NOT get initialized
	/// so call clear( defaultVal ) after reseting to initialize the grid data
	void reset( const AABB &volume, size_type resolution )
	{
		assert( resolution > 0 );

		m_volume = volume;
		m_resolution = resolution + 1;
		m_voxelSize = m_volume.getWidth() / ( m_resolution - 1 );
		m_voxels.resize( m_resolution * m_resolution * m_resolution );
	}

	void reset( const AABB &volume, mg::Real voxelSize )
	{
		const auto w = volume.getWidth();
		voxelSize = std::min( std::fabs( voxelSize ), w );

		m_volume = volume;
		m_resolution = static_cast<size_type>( w / voxelSize ) + 1;
		m_voxelSize = w / ( m_resolution - 1 );
		m_voxels.resize( m_resolution * m_resolution * m_resolution );
	}

	/// Reset values in voxels
	void clear( const value_type& defaultVal )
	{
		OMP_PARALLEL_LOOP
		for ( auto i = 0ll; i < static_cast<long long>( m_voxels.size() ); ++i )
		{
			m_voxels[i] = defaultVal;
		}
	}

	/// Insert quantity performing tri-linear interpolation
	void insertValue( const mg::Vec3D& pos, const value_type& val )
	{
		size_type i,j,k;
		mg::Vec3D localPos;
		findVoxel( pos, i, j, k, localPos );

		assert( getVoxelIdx( i+1, j+1, k+1 ) < m_voxels.size() );

		m_voxels[getVoxelIdx( i  , j  , k   )] += (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * val;
		m_voxels[getVoxelIdx( i+1, j  , k   )] += (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * val;
		m_voxels[getVoxelIdx( i  , j+1, k   )] += (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * val;
		m_voxels[getVoxelIdx( i+1, j+1, k   )] += (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * val;
		m_voxels[getVoxelIdx( i  , j  , k+1 )] +=      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * val;
		m_voxels[getVoxelIdx( i+1, j  , k+1 )] +=      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * val;
		m_voxels[getVoxelIdx( i  , j+1, k+1 )] +=      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * val;
		m_voxels[getVoxelIdx( i+1, j+1, k+1 )] +=      localPos[2]  *      localPos[1]  *      localPos[0]  * val;
	}

	/// Estimate quantity at a point performing tri-linear interpolation
	value_type estimateValueAt( const mg::Vec3D& pos ) const
	{
		size_type i,j,k;
		mg::Vec3D localPos;
		findVoxel(pos, i, j, k, localPos);

		assert( getVoxelIdx( i+1, j+1, k+1 ) < m_voxels.size() );

		return (1 - localPos[2]) * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k )]
		     + (1 - localPos[2]) * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k )]
		     + (1 - localPos[2]) *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k )]
		     + (1 - localPos[2]) *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k )]
		     +      localPos[2]  * (1 - localPos[1]) * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j  , k+1 )]
		     +      localPos[2]  * (1 - localPos[1]) *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j  , k+1 )]
		     +      localPos[2]  *      localPos[1]  * (1 - localPos[0]) * m_voxels[getVoxelIdx( i  , j+1, k+1 )]
		     +      localPos[2]  *      localPos[1]  *      localPos[0]  * m_voxels[getVoxelIdx( i+1, j+1, k+1 )];
	}

	/// Returns the voxel/sample which has min (x,y,z) coords among the nearest 8 lattice points
	/// surrounding pos. After the call o_i, o_j, o_k contain the voxel indices while o_localPos
	/// contains normalized local coords (i.e. between 0 and 1) used for tri-linear interpolation
	///
	/// TODO FIX: if pos is outside the VoxelGrid volume the nearest 8 lattice points on the grid
	/// are considered of which the min one is returned. This can lead to incorrectt results as
	/// quantities for such points being inserted accumulate on the grid's boundry
	/// or estimated for all points outside the grid are same as those on the boundry
	void findVoxel( const mg::Vec3D& pos, size_type& o_i, size_type& o_j, size_type &o_k, mg::Vec3D& o_localPos ) const
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

		o_i = static_cast<size_type>( voxelIdx[0] );
		o_j = static_cast<size_type>( voxelIdx[1] );
		o_k = static_cast<size_type>( voxelIdx[2] );
	}

private:
	constexpr inline size_type getVoxelIdx( size_type i, size_type j, size_type k ) const
	{
		return k * ( m_resolution * m_resolution ) + i * m_resolution + j;
	}

	AABB m_volume;
	mg::Real m_voxelSize;
	size_type m_resolution;
	VoxelList m_voxels;
};

using VoxelGridR = VoxelGrid<mg::Real>;
using VoxelGridVec3D = VoxelGrid<mg::Vec3D>;
