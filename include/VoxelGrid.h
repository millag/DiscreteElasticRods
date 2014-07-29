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
    unsigned insertDensity(const ngl::Vec4& pos, ngl::Real val);
//    insert particle's velovity performing trilinear interpolation
    unsigned insertVelocity(const ngl::Vec4& pos, ngl::Vec4 val);
//    calculate density at a point performing trilinear interpolation
    ngl::Real getInterpolatedDensity(const ngl::Vec4& pos) const;
//    calculate velocity at a point performing trilinear interpolation
    ngl::Vec4 getInterpolatedVelocity(const ngl::Vec4& pos) const;

protected:

    struct Voxel
    {
        ngl::Real m_density;
        ngl::Vec4 m_velocity;
    };

    const AABB& m_volume;
    unsigned m_divisions;
    unsigned m_divisionsSqr;
    ngl::Real m_voxelSize;

    std::vector<Voxel> m_voxels;

    unsigned findVoxel(const ngl::Vec4& pos, unsigned& o_i, unsigned& o_j, unsigned &o_k, ngl::Vec3& o_voxelPos) const;
    unsigned getVoxelIdx(unsigned i, unsigned j, unsigned k) const;
};

#endif // VOXELGRID_H
