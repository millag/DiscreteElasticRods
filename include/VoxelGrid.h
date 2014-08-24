#ifndef VOXELGRID_H
#define VOXELGRID_H

#include <list>
#include <vector>
#include "AABB.h"
#include "Utils.h"

class VoxelGrid {

public:
    VoxelGrid(const AABB& volume, unsigned divisions);
    ~VoxelGrid(){ }

//    allocates memory for grid
    void initialize();
//    reintialize values in cells
    void reset();
//    insert particle performing trilinear interpolation
    unsigned insertDensity(const mg::Vec3D& pos, mg::Real val);
//    insert particle's velovity performing trilinear interpolation
    unsigned insertVelocity(const mg::Vec3D& pos, const mg::Vec3D& val);
//    calculate density at a point performing trilinear interpolation
    void getInterpolatedDensity(const mg::Vec3D& pos, mg::Real& o_density) const;
//    calculate velocity at a point performing trilinear interpolation
    void getInterpolatedVelocity(const mg::Vec3D& pos, mg::Vec3D& o_velocity) const;

private:

    struct Voxel
    {
        mg::Real m_density;
        mg::Vec3D m_velocity;
    };

    const AABB& m_volume;
    unsigned m_divisions;
    unsigned m_divisionsSqr;
    mg::Real m_voxelSize;

    std::vector<Voxel> m_voxels;

    unsigned findVoxel(const mg::Vec3D& pos, unsigned& o_i, unsigned& o_j, unsigned &o_k, mg::Vec3D& o_voxelPos) const;
    unsigned getVoxelIdx(unsigned i, unsigned j, unsigned k) const;
};

#endif // VOXELGRID_H
