#include "VoxelGrid.h"

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
    typedef std::vector<Voxel>::iterator Iter;
    for (Iter it = m_voxels.begin(); it != m_voxels.end(); ++it)
    {
        it->m_density = 0;
        it->m_velocity.set(0,0,0,0);
    }
}

unsigned VoxelGrid::insertDensity(const ngl::Vec4& pos, ngl::Real val)
{
    unsigned i,j,k, nidx;
    ngl::Vec3 voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    m_voxels[ idx ].m_density += val * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i, j, k + 1);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * voxelPos.m_z;
    nidx = getVoxelIdx(i, j + 1, k);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos.m_x) * voxelPos.m_y * (1 - voxelPos.m_z);
    nidx =  getVoxelIdx(i, j + 1, k + 1);
    m_voxels[ nidx ].m_density += val * (1 - voxelPos.m_x) * voxelPos.m_y * voxelPos.m_z;

    nidx = getVoxelIdx(i + 1, j, k);
    m_voxels[ nidx ].m_density += val * voxelPos.m_x * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i + 1, j, k + 1);
    m_voxels[ nidx ].m_density += val * voxelPos.m_x * (1 - voxelPos.m_y) * voxelPos.m_z;
    nidx = getVoxelIdx(i + 1, j + 1, k);
    m_voxels[ nidx ].m_density += val * voxelPos.m_x * voxelPos.m_y * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i + 1, j + 1, k + 1);
    m_voxels[ nidx ].m_density += val * voxelPos.m_x * voxelPos.m_y * voxelPos.m_z;

    return idx;
}

unsigned VoxelGrid::insertVelocity(const ngl::Vec4& pos, ngl::Vec4 val)
{
    unsigned i, j, k, nidx;
    ngl::Vec3 voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    m_voxels[ idx ].m_velocity += val * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i, j, k + 1);;
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * voxelPos.m_z;
    nidx = getVoxelIdx(i, j + 1, k);;
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos.m_x) * voxelPos.m_y * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i, j + 1, k + 1);
    m_voxels[ nidx ].m_velocity += val * (1 - voxelPos.m_x) * voxelPos.m_y * voxelPos.m_z;
    nidx = getVoxelIdx(i + 1, j, k);
    m_voxels[ nidx ].m_velocity += val * voxelPos.m_x * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i + 1, j, k + 1);
    m_voxels[ nidx ].m_velocity += val * voxelPos.m_x * (1 - voxelPos.m_y) * voxelPos.m_z;
    nidx = getVoxelIdx(i + 1, j + 1, k);
    m_voxels[ nidx ].m_velocity += val * voxelPos.m_x * voxelPos.m_y * (1 - voxelPos.m_z);
    nidx = getVoxelIdx(i + 1, j + 1, k + 1);
    m_voxels[ nidx ].m_velocity += val * voxelPos.m_x * voxelPos.m_y * voxelPos.m_z;

    return idx;
}

ngl::Real VoxelGrid::getInterpolatedDensity(const ngl::Vec4& pos) const
{
    unsigned i,j,k;
    ngl::Vec3 voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    ngl::Real res = 0;
    res += m_voxels[ idx ].m_density * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i, j, k + 1) ].m_density * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * voxelPos.m_z;
    res += m_voxels[ getVoxelIdx(i, j + 1, k) ].m_density * (1 - voxelPos.m_x) * voxelPos.m_y * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i, j + 1, k + 1) ].m_density * (1 - voxelPos.m_x) * voxelPos.m_y * voxelPos.m_z;

    res += m_voxels[ getVoxelIdx(i + 1, j, k) ].m_density * voxelPos.m_x * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i + 1, j, k + 1) ].m_density * voxelPos.m_x * (1 - voxelPos.m_y) * voxelPos.m_z;
    res += m_voxels[ getVoxelIdx(i + 1, j + 1, k) ].m_density * voxelPos.m_x * voxelPos.m_y * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i + 1, j + 1, k + 1) ].m_density * voxelPos.m_x * voxelPos.m_y * voxelPos.m_z;

    return res;
}

ngl::Vec4 VoxelGrid::getInterpolatedVelocity(const ngl::Vec4& pos) const
{
    unsigned i,j,k;
    ngl::Vec3 voxelPos;
    unsigned idx = findVoxel(pos, i, j, k, voxelPos);
    assert(idx < m_voxels.size());

    ngl::Vec4 res(0, 0, 0, 0);
    res += m_voxels[ idx ].m_velocity * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i, j, k + 1) ].m_velocity * (1 - voxelPos.m_x) * (1 - voxelPos.m_y) * voxelPos.m_z;
    res += m_voxels[ getVoxelIdx(i, j + 1, k) ].m_velocity * (1 - voxelPos.m_x) * voxelPos.m_y * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i, j + 1, k + 1) ].m_velocity * (1 - voxelPos.m_x) * voxelPos.m_y * voxelPos.m_z;

    res += m_voxels[ getVoxelIdx(i + 1, j, k) ].m_velocity * voxelPos.m_x * (1 - voxelPos.m_y) * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i + 1, j, k + 1) ].m_velocity * voxelPos.m_x * (1 - voxelPos.m_y) * voxelPos.m_z;
    res += m_voxels[ getVoxelIdx(i + 1, j + 1, k) ].m_velocity * voxelPos.m_x * voxelPos.m_y * (1 - voxelPos.m_z);
    res += m_voxels[ getVoxelIdx(i + 1, j + 1, k + 1) ].m_velocity * voxelPos.m_x * voxelPos.m_y * voxelPos.m_z;

    return res;
}

unsigned VoxelGrid::findVoxel(const ngl::Vec4 &pos, unsigned& o_i, unsigned& o_j, unsigned &o_k, ngl::Vec3& o_voxelPos) const
{
    ngl::Vec4 origin = m_volume.getVMin();
    o_voxelPos = (pos - origin).toVec3();
    o_i = mg::clamp(std::floor(o_voxelPos.m_x / m_voxelSize), (ngl::Real)0.0, (ngl::Real)(m_divisions - 2));
    o_j = mg::clamp(std::floor(o_voxelPos.m_y / m_voxelSize), (ngl::Real)0.0, (ngl::Real)(m_divisions - 2));
    o_k = mg::clamp(std::floor(o_voxelPos.m_z / m_voxelSize), (ngl::Real)0.0, (ngl::Real)(m_divisions - 2));

    o_voxelPos.m_x = mg::clamp((ngl::Real)fmod(o_voxelPos.m_x, m_voxelSize), (ngl::Real)0, (ngl::Real)1);
    o_voxelPos.m_y = mg::clamp((ngl::Real)fmod(o_voxelPos.m_y, m_voxelSize), (ngl::Real)0, (ngl::Real)1);
    o_voxelPos.m_z = mg::clamp((ngl::Real)fmod(o_voxelPos.m_z, m_voxelSize), (ngl::Real)0, (ngl::Real)1);

    return getVoxelIdx(o_i, o_j, o_k);
}

unsigned VoxelGrid::getVoxelIdx(unsigned i, unsigned j, unsigned k) const
{
    return k * m_divisionsSqr + i * m_divisions + j;
}
