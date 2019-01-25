#pragma once

#include "AABB.h"
#include "Utils.h"

/// Regularly spaced grid encompassing some volume in 3D space. Each voxel
/// represents only a single sample point on the grid and not a volume. How
/// many voxels/samples the grid has depends on the grid resolution - f.e. a
/// VoxelGrid with non-zero volume and resolution 1 has 8 samples, one per
/// each vertex of the volume bbox. The space between each voxel is not represented
/// in the dataset, but a value for each point in 3D space is approximated via
/// tri-linear interpolation. When inserting a value in the VoxelGrid it is distributed
/// among the 8 surrounding voxels using again tri-linear interpolation
class VoxelGrid
{
public:

	const AABB& getVolume() const { return m_volume; }
	unsigned getResolution() const { return m_resolution - 1; }
	mg::Real getVolexSize() const { return m_voxelSize; }

///	Setsup and allocates memory for grid with the given volume and resolution
/// Voxel values are kept on the
	void reset( const AABB &volume, unsigned resolution );
	void reset( const AABB &volume, mg::Real voxelSize );

///	Resets values in voxels
	void clear();

///	Inserts mass particle performing tri-linear interpolation
	unsigned insertDensity( const mg::Vec3D& pos, mg::Real val );

///	Inserts mass particle velovity performing tri-linear interpolation
	unsigned insertVelocity( const mg::Vec3D& pos, const mg::Vec3D& val );

///	Returns density at a point performing tri-linear interpolation
	void getInterpolatedDensity( const mg::Vec3D& pos, mg::Real& o_density ) const;

///	Returns velocity at a point performing tri-linear interpolation
	void getInterpolatedVelocity( const mg::Vec3D& pos, mg::Vec3D& o_velocity ) const;

private:
	constexpr inline unsigned getVoxelIdx( unsigned i, unsigned j, unsigned k ) const
	{
		return k * ( m_resolution * m_resolution ) + i * m_resolution + j;
	}

	unsigned findVoxel( const mg::Vec3D& pos, unsigned& o_i, unsigned& o_j, unsigned &o_k, mg::Vec3D& o_voxelPos ) const;

	struct Voxel
	{
		mg::Real m_density;
		mg::Vec3D m_velocity;
	};

	AABB m_volume;
	unsigned m_resolution;
	mg::Real m_voxelSize;
	std::vector<Voxel> m_voxels;
};
