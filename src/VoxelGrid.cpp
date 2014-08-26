#include "VoxelGrid.h"
#include <cassert>

VoxelGrid::VoxelGrid(const AABB &volume, unsigned divisions):
    m_volume(volume), m_divisions(divisions + 1), m_divisionsSqr((divisions + 1) * (divisions + 1))
{
    assert(m_divisions > 1);
    m_voxelSize = m_volume.getWidth() / (m_divisions - 1);
}

void VoxelGrid::initialize()
{
    m_voxels.resize(m_divisionsSqr * m_divisions);
    reset();
}

void VoxelGrid::reset()
{
    m_voxelSize = m_volume.getWidth() / (m_divisions - 1);

#ifdef MULTI_THREAD
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < m_voxels.size(); ++i)
    {
        m_voxels[i].m_density = 0;
        m_voxels[i].m_velocity.zero();
    }
}

unsigned VoxelGrid::insertDensity(const mg::Vec3D& pos, mg::Real val)
{
    unsigned i,j,k, nidx;
    mg::Vec3D voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    m_voxels[ idx ].m_density += val * (1 - voxelPos[0]) * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i, j, k + 1);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos[0]) * (1 - voxelPos[1]) * voxelPos[2];
    nidx = getVoxelIdx(i, j + 1, k);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos[0]) * voxelPos[1] * (1 - voxelPos[2]);
    nidx =  getVoxelIdx(i, j + 1, k + 1);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos[0]) * voxelPos[1] * voxelPos[2];

    nidx = getVoxelIdx(i + 1, j, k);
    m_voxels[ nidx ].m_density += val * voxelPos[0] * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i + 1, j, k + 1);
    m_voxels[ nidx ].m_density += val * voxelPos[0] * (1 - voxelPos[1]) * voxelPos[2];
    nidx = getVoxelIdx(i + 1, j + 1, k);
    m_voxels[ nidx ].m_density += val * voxelPos[0] * voxelPos[1] * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i + 1, j + 1, k + 1);
    m_voxels[ nidx ].m_density += val * voxelPos[0] * voxelPos[1] * voxelPos[2];

    return idx;
}

unsigned VoxelGrid::insertVelocity(const mg::Vec3D& pos, const mg::Vec3D& val)
{
    unsigned i, j, k, nidx;
    mg::Vec3D voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    m_voxels[ idx ].m_velocity += val * (1 - voxelPos[0]) * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i, j, k + 1);;
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos[0]) * (1 - voxelPos[1]) * voxelPos[2];
    nidx = getVoxelIdx(i, j + 1, k);;
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos[0]) * voxelPos[1] * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i, j + 1, k + 1);
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos[0]) * voxelPos[1] * voxelPos[2];
    nidx = getVoxelIdx(i + 1, j, k);
    m_voxels[ nidx ].m_velocity += val * voxelPos[0] * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i + 1, j, k + 1);
    m_voxels[ nidx ].m_velocity += val * voxelPos[0] * (1 - voxelPos[1]) * voxelPos[2];
    nidx = getVoxelIdx(i + 1, j + 1, k);
    m_voxels[ nidx ].m_velocity += val * voxelPos[0] * voxelPos[1] * (1 - voxelPos[2]);
    nidx = getVoxelIdx(i + 1, j + 1, k + 1);
    m_voxels[ nidx ].m_velocity += val * voxelPos[0] * voxelPos[1] * voxelPos[2];

    return idx;
}

void VoxelGrid::getInterpolatedDensity(const mg::Vec3D& pos, mg::Real& o_density) const
{
    unsigned i,j,k;
    mg::Vec3D voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    o_density = 0;
    o_density += m_voxels[ idx ].m_density * (1 - voxelPos[0]) * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    o_density += m_voxels[ getVoxelIdx(i, j, k + 1) ].m_density * (1 - voxelPos[0]) * (1 - voxelPos[1]) * voxelPos[2];
    o_density += m_voxels[ getVoxelIdx(i, j + 1, k) ].m_density * (1 - voxelPos[0]) * voxelPos[1] * (1 - voxelPos[2]);
    o_density += m_voxels[ getVoxelIdx(i, j + 1, k + 1) ].m_density * (1 - voxelPos[0]) * voxelPos[1] * voxelPos[2];

    o_density += m_voxels[ getVoxelIdx(i + 1, j, k) ].m_density * voxelPos[0] * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    o_density += m_voxels[ getVoxelIdx(i + 1, j, k + 1) ].m_density * voxelPos[0] * (1 - voxelPos[1]) * voxelPos[2];
    o_density += m_voxels[ getVoxelIdx(i + 1, j + 1, k) ].m_density * voxelPos[0] * voxelPos[1] * (1 - voxelPos[2]);
    o_density += m_voxels[ getVoxelIdx(i + 1, j + 1, k + 1) ].m_density * voxelPos[0] * voxelPos[1] * voxelPos[2];
}

void VoxelGrid::getInterpolatedVelocity(const mg::Vec3D& pos, mg::Vec3D& o_velocity) const
{
    unsigned i,j,k;
    mg::Vec3D voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    o_velocity.zero();
    o_velocity += m_voxels[ idx ].m_velocity * (1 - voxelPos[0]) * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    o_velocity += m_voxels[ getVoxelIdx(i, j, k + 1) ].m_velocity * (1 - voxelPos[0]) * (1 - voxelPos[1]) * voxelPos[2];
    o_velocity += m_voxels[ getVoxelIdx(i, j + 1, k) ].m_velocity * (1 - voxelPos[0]) * voxelPos[1] * (1 - voxelPos[2]);
    o_velocity += m_voxels[ getVoxelIdx(i, j + 1, k + 1) ].m_velocity * (1 - voxelPos[0]) * voxelPos[1] * voxelPos[2];

    o_velocity += m_voxels[ getVoxelIdx(i + 1, j, k) ].m_velocity * voxelPos[0] * (1 - voxelPos[1]) * (1 - voxelPos[2]);
    o_velocity += m_voxels[ getVoxelIdx(i + 1, j, k + 1) ].m_velocity * voxelPos[0] * (1 - voxelPos[1]) * voxelPos[2];
    o_velocity += m_voxels[ getVoxelIdx(i + 1, j + 1, k) ].m_velocity * voxelPos[0] * voxelPos[1] * (1 - voxelPos[2]);
    o_velocity += m_voxels[ getVoxelIdx(i + 1, j + 1, k + 1) ].m_velocity * voxelPos[0] * voxelPos[1] * voxelPos[2];
}

unsigned VoxelGrid::findVoxel(const mg::Vec3D& pos, unsigned& o_i, unsigned& o_j, unsigned &o_k, mg::Vec3D& o_voxelPos) const
{
    o_voxelPos = pos - m_volume.getVMin();
    o_i = mg::clamp(std::floor(o_voxelPos[0] / m_voxelSize), (mg::Real)0.0, (mg::Real)(m_divisions - 2));
    o_j = mg::clamp(std::floor(o_voxelPos[1] / m_voxelSize), (mg::Real)0.0, (mg::Real)(m_divisions - 2));
    o_k = mg::clamp(std::floor(o_voxelPos[2] / m_voxelSize), (mg::Real)0.0, (mg::Real)(m_divisions - 2));

    o_voxelPos[0] = mg::clamp((mg::Real)fmod(o_voxelPos[0], m_voxelSize), (mg::Real)0, (mg::Real)1);
    o_voxelPos[1] = mg::clamp((mg::Real)fmod(o_voxelPos[1], m_voxelSize), (mg::Real)0, (mg::Real)1);
    o_voxelPos[2] = mg::clamp((mg::Real)fmod(o_voxelPos[2], m_voxelSize), (mg::Real)0, (mg::Real)1);

    return getVoxelIdx(o_i, o_j, o_k);
}

unsigned VoxelGrid::getVoxelIdx(unsigned i, unsigned j, unsigned k) const
{
    return k * m_divisionsSqr + i * m_divisions + j;
}
