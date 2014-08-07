#include "ElasticRod.h"

#include <cmath>

#include "Utils.h"


ElasticRod::ElasticRod()
{
    m_bendStiffnes = 0.9;
    m_twistStiffness = 1.0;
    m_nIterations = 4;
    m_maxElasticForceThreshold = 100;

    m_B.identity();
    m_B *= m_bendStiffnes;
}

ElasticRod::~ElasticRod()
{ }


void ElasticRod::init(const std::vector<mg::Vec3D>& restpos,
          const mg::Vec3D u0,
          const std::vector<mg::Vec3D>& pos,
          const std::vector<mg::Vec3D>& vel,
          const std::vector<mg::Real>& mass,
          const std::vector<mg::Real>& twistAngle,
          const std::vector<bool>& isFixed)
{
    assert( pos.size() > 2 );
    assert( pos.size() == restpos.size() );
    assert( pos.size() == vel.size() );
    assert( pos.size() == mass.size() );
    assert( pos.size() == isFixed.size() );
    assert( twistAngle.size() == (pos.size() - 1) );

    m_totalTwist = 8.;

    m_u0 = u0;
    m_ppos = pos;
    m_pvel = vel;
    m_pmass = mass;
    m_twistAngle = twistAngle;
    m_pIsFixed = isFixed;

    m_edges.resize(restpos.size() - 1);

    m_restEdgeL.resize(m_edges.size());
    m_restRegionL.resize(m_edges.size());
    m_restWprev.resize(m_edges.size());
    m_restWnext.resize(m_edges.size());

    m_kb.resize(m_edges.size());
    m_m1.resize(m_edges.size());
    m_m2.resize(m_edges.size());

//    compute edges & lengths for rest shape
    computeEdges(restpos, m_edges);
    computeLengths(m_edges, m_restEdgeL, m_restRegionL);
//    compute kb and material frame (Bishop frame) for rest shape
    computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);
//    precompute material curvature for rest shape
    computeMaterialCurvature(m_kb, m_m1, m_m2, m_restWprev,  m_restWnext);

//    m_u = m_m1;

//    initialize current state
    updateCurrentState();

//    computeEdges(m_ppos, m_edges);
//    computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);
//    computeMaterialFrame(m_twistAngle, m_m1, m_m2);
}


void ElasticRod::computeEdges(const std::vector<mg::Vec3D>& vertices,
                              std::vector<mg::Vec3D>& o_edges) const
{
    for (unsigned i = 0; i < vertices.size() - 1; ++i)
    {
        o_edges[i] = vertices[i + 1] - vertices[i];
        assert( std::isfinite(mg::length_squared(o_edges[i])) );
    }
}

void ElasticRod::computeLengths(const std::vector<mg::Vec3D>& edges,
                                std::vector<mg::Real> &o_edgeL,
                                std::vector<mg::Real> &o_regionL) const
{
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        o_edgeL[i] = edges[i].length();
        assert( std::isfinite(o_edgeL[i]) && o_edgeL[i] > mg::ERR );
    }

    o_regionL[0] = 0.0;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_regionL[i] = o_edgeL[i - 1] + o_edgeL[i];
    }
}

void ElasticRod::computeMaterialCurvature(const std::vector<mg::Vec3D>& kb,
                                        const std::vector<mg::Vec3D>& m1,
                                        const std::vector<mg::Vec3D>& m2,
                                        std::vector<mg::Vec2D>& o_Wprev,
                                        std::vector<mg::Vec2D>& o_Wnext) const
{
    o_Wprev[0].zero();
    o_Wnext[0].zero();
    for (unsigned i = 1; i < kb.size(); ++i)
    {
        computeW(kb[i], m1[i - 1], m2[i - 1], o_Wprev[i]);
        computeW(kb[i], m1[i], m2[i], o_Wnext[i]);
    }
}

void ElasticRod::computeW(const mg::Vec3D& kb, const mg::Vec3D& m1, const mg::Vec3D& m2, mg::Vec2D& o_wij) const
{
    o_wij.set( mg::dot( kb, m2 ), -mg::dot( kb, m1 ));
}

void ElasticRod::computeKB(const std::vector<mg::Vec3D>& edges,
                           std::vector<mg::Vec3D>& o_kb) const
{
    o_kb[0].zero();

    double magnitude;
    double sinPhi, cosPhi;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
//     NOTE: as metioned in the paper the following formula produces |kb| = 2 * tan(phi / 2),
//     where phi is the angle of rotation defined by the 2 edges
        o_kb[i] = 2 * mg::cross( edges[i - 1], edges[i] ) /
                (m_restEdgeL[i - 1] * m_restEdgeL[i] + mg::dot( edges[i - 1], edges[i]) );

//      following for testing purpose only - need to assert that |kb| = 2 * tan(phi / 2)
        magnitude = mg::length_squared(o_kb[i]);
        extractSinAndCos(magnitude, sinPhi, cosPhi);
    }
}

void ElasticRod::extractSinAndCos(const double &magnitude,
                                  double &o_sinPhi, double &o_cosPhi) const
{

    o_cosPhi = mg::sqrt_safe(4.0 / (4.0 + magnitude));
    o_sinPhi = mg::sqrt_safe(magnitude / (4.0 + magnitude));

//    actually tan maps [-pi/2, pi/2] to [-inf, inf] so asserts don't check correcly
//    interestingly for arbitrary real number a bounded angle can be extracted
    assert( std::isfinite(magnitude) && magnitude >= 0.0  && std::isfinite(o_sinPhi) && std::isfinite(o_cosPhi) );
    assert( fabs((o_sinPhi * o_sinPhi + o_cosPhi * o_cosPhi) - 1.0) < mg::ERR );
    assert( fabs(magnitude - 4.0 * o_sinPhi * o_sinPhi / (o_cosPhi * o_cosPhi)) < mg::ERR );
}

void ElasticRod::computeBishopFrame(const mg::Vec3D& u0,
                        const std::vector<mg::Vec3D>& edges,
                        std::vector<mg::Vec3D>& o_kb,
                        std::vector<mg::Vec3D>& o_u,
                        std::vector<mg::Vec3D>& o_v) const
{
//    check if Bishop frame is correctly defined
    assert( fabs(mg::length_squared( u0 ) - 1) < mg::ERR );
//    std::cout << "ANGLE u0 edge[0] " << mg::deg( std::acos( mg::dot( u0, mg::normalize(edges[0]) ) ) ) << std::endl;

    computeKB(edges, o_kb);

    o_u[0] = u0;
    o_v[0] = mg::cross(edges[0], o_u[0]);
    o_v[0].normalize();

    double magnitude;
    double sinPhi, cosPhi;
//    compute Bishop frame for current configuration given u0 and edges
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        magnitude = mg::length_squared(o_kb[i]);
        if (magnitude < mg::ERR )
        {
            o_u[i] = o_u[i - 1];
            o_v[i] = o_v[i - 1];
            continue;
        }

//        here sinPhi and cosPhi are derived from the length of kb |kb| = 2 * tan( phi/2 )
        extractSinAndCos(magnitude, sinPhi, cosPhi);

//        rotate frame u axis around kb
        mg::Quaternion q(cosPhi, sinPhi * mg::normalize(o_kb[i]));
        mg::Quaternion p(0, o_u[i - 1]);
        p = q * p * mg::conjugate(q);

        o_u[i].set(p[1], p[2], p[3]);
        o_u[i].normalize();
        o_v[i] = mg::cross(edges[i], o_u[i]);
        o_v[i].normalize();

        assert( mg::length_squared(o_kb[i]) > mg::ERR );
        assert( fabs(p[0]) <  mg::ERR );
        assert( fabs(mg::length_squared(o_u[i]) - 1) < mg::ERR );
        assert( fabs(mg::length_squared(o_v[i]) - 1) < mg::ERR );
//        std::cout << "LENGTH u[" << i << "] " << mg::length_squared(o_u[i]) << std::endl;
//        std::cout << "ANGLE u edge[" << i << "] " << mg::deg( std::acos( mg::dot( o_u[i], mg::normalize(edges[i]) ) ) ) << std::endl;
//        std::cout << "ANGLE v edge[" << i << "] " << mg::deg( std::acos( mg::dot( o_v[i], mg::normalize(edges[i]) ) ) ) << std::endl;
    }
}

void ElasticRod::computeMaterialFrame(const std::vector<mg::Real>& twistAngle,
                                      std::vector<mg::Vec3D>& io_m1,
                                      std::vector<mg::Vec3D>& io_m2) const
{
    mg::Real sinQ, cosQ;
    mg::Vec3D m1, m2;
    for (unsigned i = 0; i < io_m1.size(); ++i)
    {
        cosQ = std::cos(twistAngle[i]);
        sinQ = std::sqrt(1 - cosQ * cosQ);

        m1 = cosQ * io_m1[i] + sinQ * io_m2[i];
        m2 = -sinQ * io_m1[i] + cosQ * io_m2[i];

        io_m1[i] = m1;
        io_m2[i] = m2;

        assert( fabs( mg::length_squared( m1 ) - 1 ) < mg::ERR );
        assert( fabs( mg::length_squared( m2 ) - 1 ) < mg::ERR);
        assert( fabs( mg::dot( m1, m2 )) < mg::ERR );
//        std::cout << "ANGLE m1 edge[" << i << "] " << mg::deg( std::acos( mg::dot( m1, mg::normalize(edges[i]) ) ) ) << std::endl;
//        std::cout << "ANGLE m2 edge[" << i << "] " << mg::deg( std::acos( mg::dot( m2, mg::normalize(edges[i]) ) ) ) << std::endl;
//        assert( fabs( mg::dot( m1, mg::normalize(edges[i]) )) < mg::THRESHODL );
//        assert( fabs( mg::dot( m2, mg::normalize(edges[i]) )) < mg::THRESHODL );
    }
}

/* Use Verlet integration. */
void ElasticRod::update(mg::Real dt)
{
//    integrate centerline
    std::vector<mg::Vec3D> forces(m_ppos.size(), mg::Vec3D(0,0,0));
    computeForces(m_ppos, forces);

    std::vector<mg::Vec3D> prevPos(m_ppos.size());
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        prevPos[i] = m_ppos[i];
        m_pvel[i] += dt * forces[i] / m_pmass[i];
        m_ppos[i] += m_pvel[i] * dt;
    }

//    solve constraints
    mg::Vec3D e;
    mg::Real l, l1, l2;
    for (unsigned k = 0; k < m_nIterations; ++k)
    {
        for (unsigned i = 0; i < m_ppos.size() - 1; ++i)
        {
            if (m_pIsFixed[i] && m_pIsFixed[i + 1])
            {
                continue;
            }

            e = m_ppos[i + 1] - m_ppos[i];
            l = e.length();
            if (l > mg::ERR)
            {
                e.normalize();
            }
            l = (l - m_restEdgeL[i]);

            if (m_pIsFixed[i])
            {
                l1 = 0;
                l2 = -l;
            }
            else if (m_pIsFixed[i + 1])
            {
                l1 = l;
                l2 = 0;
            }
            else
            {
                l1 = m_pmass[i + 1] / (m_pmass[i] + m_pmass[i + 1]) * l;
                l2 = -m_pmass[i] / (m_pmass[i] + m_pmass[i + 1]) * l;
            }

            m_ppos[i] += l1 * e;
            m_ppos[i + 1] += l2 * e;
        }
    }

//    update velocities:
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        m_pvel[i] = (m_ppos[i] - prevPos[i]) / dt;
    }

    updateCurrentState();
}

void ElasticRod::computeForces(const std::vector<mg::Vec3D>& vertices, std::vector<mg::Vec3D>& o_forces)
{
    std::vector<mg::Vec3D> externalForces(vertices.size(), mg::Vec3D(0,0,0));
    std::vector<mg::Vec3D> elasticForces(vertices.size(), mg::Vec3D(0,0,0));

    computeExternalForces(vertices, externalForces);
    computeElasticForces(vertices, elasticForces);

    for (unsigned i = 0; i < o_forces.size(); ++i)
    {
        if (m_pIsFixed[i])
        {
            continue;
        }

        o_forces[i] += externalForces[i] + elasticForces[i];
    }
}


void ElasticRod::computeExternalForces(const std::vector<mg::Vec3D>& vertices,
                                       std::vector<mg::Vec3D>& o_forces) const
{
    for (unsigned i = 0; i < vertices.size(); ++i)
    {
        if (m_pIsFixed[i])
        {
            continue;
        }
        o_forces[i] += mg::GRAVITY * m_pmass[i];
    }
}

void ElasticRod::computeElasticForces(const std::vector<mg::Vec3D>& vertices,
                                      std::vector<mg::Vec3D>& o_forces)
{
//    TODO: need to implement minimization of energy for twistAngles first
    std::vector<mg::Matrix3D> minusGKB(m_edges.size());
    std::vector<mg::Matrix3D> plusGKB(m_edges.size());
    std::vector<mg::Matrix3D> eqGKB(m_edges.size());
    computeGradientKB(m_kb, m_edges, minusGKB, plusGKB, eqGKB);

    std::vector<mg::Vec3D> minusGH(m_kb.size());
    std::vector<mg::Vec3D> plusGH(m_kb.size());
    std::vector<mg::Vec3D> eqGH(m_kb.size());
    computeGradientHolonomy(m_kb, minusGH, plusGH, eqGH);

    mg::Vec2D wkj;
    mg::Matrix2D J;
    mg::matrix_rotation_2D(J, mg::Constants::pi_over_2());

//    compute dE/d(theta_n)
    unsigned n = m_edges.size() - 1;
    mg::Real dEdQn;
    computedEdQj(n, J * m_B, dEdQn);

    mg::Matrix23D GW;
    mg::Vec3D GH;
    for (unsigned i = 0; i < vertices.size(); ++i)
    {
        if (m_pIsFixed[i])
        {
            continue;
        }

        for (unsigned k = std::max((int)(i) - 1, 1); k < m_edges.size(); ++k)
        {
            computeW(m_kb[k], m_m1[k - 1], m_m2[k - 1], wkj);
            computeGradientCurvature(i, k, k - 1,
                                     minusGKB, plusGKB, eqGKB,
                                     minusGH, plusGH, eqGH,
                                     wkj,
                                     J,
                                     GW);
            o_forces[i] -= mg::transpose(GW) * (m_B * (wkj - m_restWprev[k]) / m_restRegionL[k]);
            assert( std::isfinite(mg::length_squared(o_forces[i])) );

            computeW(m_kb[k], m_m1[k], m_m2[k], wkj);
            computeGradientCurvature(i, k, k,
                                     minusGKB, plusGKB, eqGKB,
                                     minusGH, plusGH, eqGH,
                                     wkj,
                                     J,
                                     GW);
            o_forces[i] -= mg::transpose(GW) * (m_B * (wkj - m_restWnext[k]) / m_restRegionL[k]);
            assert( std::isfinite(mg::length_squared(o_forces[i])) );
        }

////        add  dE/d(theta_n) * gradient holonomy(GH)
        computeGradientHolonomySum(i, n, minusGH, plusGH, eqGH, GH);
        o_forces[i] += dEdQn * GH;

        assert( std::isfinite(mg::length_squared(o_forces[i])) );

//        need to limit the force otherwise when e[i] ~ -e[i - 1] | kb | goes to infinity => force goes to infinity
        if (mg::length_squared(o_forces[i]) > m_maxElasticForceThreshold * m_maxElasticForceThreshold)
        {
            o_forces[i].normalize();
            o_forces[i] *= m_maxElasticForceThreshold;
        }
    }
}

//    mg::Real dEdQj;
//        for (unsigned j = std::max((int)(i) - 1, 0); j < m_edges.size(); ++j)
//        {
//            computedEdQj(j, J * m_B, dEdQj);
//            computeGradientHolonomySum(i, n, minusGH, plusGH, eqGH, GH);
//            o_forces[i] += dEdQj * GH;
//        }


void ElasticRod::computeGradientKB(const std::vector<mg::Vec3D> &kb,
                                   const std::vector<mg::Vec3D> &edges,
                                   std::vector<mg::Matrix3D>& o_minusGKB,
                                   std::vector<mg::Matrix3D>& o_plusGKB,
                                   std::vector<mg::Matrix3D>& o_eqGKB) const
{
    std::vector<mg::Matrix3D> edgeMatrix(edges.size());
    computeEdgeMatrices(edges, edgeMatrix);

    o_minusGKB[0].zero();
    o_plusGKB[0].zero();
    o_eqGKB[0].zero();

    mg::Real scalarFactor;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        scalarFactor = (m_restEdgeL[i - 1] * m_restEdgeL[i] + mg::dot( edges[i - 1], edges[i] ));
        assert( std::isfinite(scalarFactor) && fabs(scalarFactor) > 0 );

        o_minusGKB[i] = (2.0 * edgeMatrix[i] + mg::outer( kb[i], edges[i] )) / scalarFactor;
        o_plusGKB[i] = (2.0 * edgeMatrix[i - 1] - mg::outer( kb[i], edges[i - 1] )) / scalarFactor;
        o_eqGKB[i] = -(o_minusGKB[i] + o_plusGKB[i]);
    }
}

void ElasticRod::computeEdgeMatrices(const std::vector<mg::Vec3D>& edges,
                                     std::vector<mg::Matrix3D>& o_edgeMatrices) const
{
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        mg::matrix_skew_symmetric(o_edgeMatrices[i], edges[i]);
    }
}

void ElasticRod::computeGradientHolonomy(const std::vector<mg::Vec3D> &kb,
                                       std::vector<mg::Vec3D>& o_minusGH,
                                       std::vector<mg::Vec3D>& o_plusGH,
                                       std::vector<mg::Vec3D>& o_eqGH) const
{
    o_minusGH[0].zero();
    o_plusGH[0].zero();
    o_eqGH[0].zero();

    for (unsigned i = 1; i < kb.size(); ++i)
    {
        o_minusGH[i] = kb[i] * 0.5 / m_restEdgeL[i - 1];
        o_plusGH[i]  = kb[i] * -0.5 / m_restEdgeL[i];
        o_eqGH[i] = -(o_minusGH[i] + o_plusGH[i]);

        assert( std::isfinite( mg::length_squared(o_minusGH[i]) ) );
        assert( std::isfinite( mg::length_squared(o_plusGH[i]) ) );
        assert( std::isfinite( mg::length_squared(o_eqGH[i]) ) );
    }
}

void ElasticRod::computeGradientCurvature(unsigned i, unsigned k, unsigned j,
                                        const std::vector<mg::Matrix3D>& minusGKB,
                                        const std::vector<mg::Matrix3D>& plusGKB,
                                        const std::vector<mg::Matrix3D>& eqGKB,
                                        const std::vector<mg::Vec3D>& minusGH,
                                        const std::vector<mg::Vec3D>& plusGH,
                                        const std::vector<mg::Vec3D>& eqGH,
                                        const mg::Vec2D &wkj,
                                        const mg::Matrix2D J,
                                        mg::Matrix23D &o_GW) const
{
    assert(k >= (i - 1) && (j == k || j == (k - 1)) && j < m_m1.size());

//    need to make o_GW zero 3x2 matrix
    o_GW.zero();
//    compute gradient KB(GKB) term
    if (k < i + 2)
    {
        o_GW(0,0) = m_m2[j][0];
        o_GW(0,1) = m_m2[j][1];
        o_GW(0,2) = m_m2[j][2];

        o_GW(1,0) = -m_m1[j][0];
        o_GW(1,1) = -m_m1[j][1];
        o_GW(1,2) = -m_m1[j][2];

        if (k == i - 1)
        {
            o_GW = o_GW * plusGKB[k];
        } else if (k == i)
        {
            o_GW = o_GW * eqGKB[k];
        } else if (k == i + 1)
        {
            o_GW = o_GW * minusGKB[k];
        }
    }
//    compute gradient Holonomy(GH) term
    mg::Vec3D GH;
    computeGradientHolonomySum(i, j, minusGH, plusGH, eqGH, GH);
    o_GW -= J * mg::outer(wkj, GH);
}

void ElasticRod::computeGradientHolonomySum(unsigned i , unsigned j,
                                 const std::vector<mg::Vec3D>& minusGH,
                                 const std::vector<mg::Vec3D>& plusGH,
                                 const std::vector<mg::Vec3D>& eqGH,
                                 mg::Vec3D& o_GH) const
{
    o_GH.zero();

    if (j >= (i - 1) && i > 1 && (i - 1) < plusGH.size())
    {
        o_GH += plusGH[i - 1];
    }
    if (j >= i && i < (eqGH.size() - 1))
    {
        o_GH += eqGH[i];
    }
    if (j >= (i + 1) && i < (minusGH.size() - 2))
    {
        o_GH += minusGH[i + 1];
    }
}

void ElasticRod::parallelTransportFrame(const mg::Vec3D& e0, const mg::Vec3D& e1,
                                        mg::Vec3D& io_u) const
{
    mg::Vec3D axis = 2 * mg::cross(e0, e1) /
            (m_restEdgeL[0] * m_restEdgeL[0] + mg::dot(e0, e1));

    double magnitude = mg::length_squared(axis);
    if (magnitude < mg::ERR)
    {
        axis = mg::cross(e1, io_u);
        io_u = mg::cross(axis, e1);
        io_u.normalize();
        return;
    }

    double sinPhi, cosPhi;
    extractSinAndCos(magnitude, sinPhi, cosPhi);

    mg::Quaternion q(cosPhi, sinPhi * mg::normalize(axis));
    mg::Quaternion p(0, io_u);
    p = q * p * mg::conjugate(q);

    io_u.set(p[1], p[2], p[3]);
    io_u.normalize();
    assert(fabs(mg::length_squared(io_u) - 1) < mg::ERR);
}

void ElasticRod::updateCurrentState()
{

//    parallel transport first frame in time
    mg::Vec3D e0 = m_edges[0];
    computeEdges(m_ppos, m_edges);
    mg::Vec3D e1 = m_edges[0];
    parallelTransportFrame(e0, e1, m_u0);

//    std::cout << "ANGLE u0 edge[0] " << mg::deg( std::acos( mg::dot(m_u0, t1) ) ) << std::endl;

    computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);
    computeMaterialFrame(m_twistAngle, m_m1, m_m2);
}


void ElasticRod::computedEdQj(unsigned j, const mg::Matrix2D JB, mg::Real& o_dEQj) const
{
    assert( j < m_kb.size() );

    o_dEQj = 0.0;
//    compute first term dWj/dQj + 2 * beta * mj / lj
    mg::Vec2D wij;
    if (j > 0)
    {
        computeW(m_kb[j], m_m1[j], m_m2[j], wij);
        o_dEQj += mg::dot(wij, JB * (wij - m_restWnext[j]));
        o_dEQj += 2 * m_twistStiffness * (m_twistAngle[j] - m_twistAngle[j - 1]);
        o_dEQj /= m_restRegionL[j];
    }
//    compute second term dWj+1/dQj + 2 * beta * mj+1 / lj+1
    if (j < m_edges.size() - 1)
    {
        computeW(m_kb[j + 1], m_m1[j], m_m2[j], wij);
        o_dEQj += mg::dot(wij, JB * (wij - m_restWprev[j + 1]));
        o_dEQj += 2 * m_twistStiffness * (m_twistAngle[j + 1] - m_twistAngle[j]);
        o_dEQj /= m_restRegionL[j + 1];
    }
}

//============================  =======================


//void ElasticRod::updateCurrentState()
//{
//    std::vector<mg::Vec3D> edges_prev = m_edges;
//    computeEdges(m_ppos, m_edges);


////    parallel transport Bishop frame in time
//    for (unsigned i = 0; i < m_edges.size(); ++i)
//    {
//        parallelTransportFrame(edges_prev[i], m_edges[i], m_u[i]);
//    }
//    m_u0 = m_u[0];

////    std::cout << "ANGLE u0 edge[0] " << mg::deg( std::acos( mg::dot(m_u0, t1) ) ) << std::endl;
//    computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);

//    mg::Real cosQ, sinQ;
//    m_twistAngle[0] = 0;
//    for (unsigned i = 1; i < m_edges.size(); ++i)
//    {
//        cosQ = mg::dot(m_u[i], m_m1[i]);
//        sinQ = mg::dot(m_u[i], m_m2[i]);
//        m_twistAngle[i] = mg::sign(sinQ) * mg::acos_safe(cosQ);
////        std::cout<< "Twist angle [" << i << "] = " << cosQ * cosQ + sinQ * sinQ << std::endl;
////        std::cout<< "Twist angle [" << i << "] = " << mg::deg(m_twistAngle[i]) << std::endl;
//    }

//    computeMaterialFrame(m_twistAngle, m_m1, m_m2);

//    for (unsigned i = 1; i < m_edges.size(); ++i)
//    {
//        std::cout<< "Twist angle [" << i << "] = " << mg::dot( m_u[i], m_m1[i] ) << std::endl;
//    }
//    m_u = m_m1;
//}




//void ElasticRod::computeElasticForces(const std::vector<mg::Vec3D>& vertices, std::vector<mg::Vec3D>& o_forces)
//{
////    TODO: need to implement minimization of energy for twistAngles first

//    std::vector<mg::Matrix3D> minusGKB(3);
//    std::vector<mg::Matrix3D> eqGKB(3);
//    std::vector<mg::Matrix3D> plusGKB(3);
//    std::vector<mg::Vec3D> minusGH(3);
//    std::vector<mg::Vec3D> eqGH(3);
//    std::vector<mg::Vec3D> plusGH(3);
//    for (unsigned i = 0; i < minusGKB.size(); ++i)
//    {
//        minusGKB.zero();
//        eqGKB.zero();
//        plusGKB.zero();
//        minusGH.zero();
//        eqGH.zero();
//        plusGH.zero();
//    }



//    mg::Matrix23D gw;
//    for (unsigned i = 0; i < vertices.size(); ++i)
//    {
//        computeGKBandGH(i + 1, minusGKB, eqGKB, plusGKB, minusGH, eqGH, plusGH);

//        if (m_pIsFixed[i])
//        {
//            continue;
//        } ...



//void ElasticRod::computeGKBandGH(unsigned idx,
//                                std::vector<mg::Matrix3D> &o_minusGKB,
//                                std::vector<mg::Matrix3D> &o_eqGKB,
//                                std::vector<mg::Matrix3D> &o_plusGKB,
//                                std::vector<mg::Vec3D> &o_minusGH,
//                                std::vector<mg::Vec3D> &o_eqGH,
//                                std::vector<mg::Vec3D> &o_plusGH) const
//{
//    assert(idx > 0 && idx < m_edges.size());

//    for (unsigned i = 0; i < o_minusGKB.size() - 1; ++i)
//    {
//        o_minusGKB[i] = o_minusGKB[i + 1];
//        o_eqGKB[i] = o_eqGKB[i + 1];
//        o_plusGKB[i] = o_plusGKB[i + 1];

//        o_minusGH[i] = o_minusGH[i + 1];
//        o_eqGH[i] = o_eqGH[i + 1];
//        o_plusGH[i] = o_plusGH[i + 1];
//    }

//    mg::Real scalarFactor = (m_restEdgeL[idx - 1] * m_restEdgeL[idx] + mg::dot( m_edges[idx - 1], m_edges[idx] ));;
//    assert( std::isfinite(scalarFactor) && fabs(scalarFactor) > 0 );
//    mg::Matrix3D edgeMatrix;
//    mg::matrix_skew_symmetric(edgeMatrix, edges[idx]);

//    o_minusGKB[i] = (2.0 * edgeMatrix[i] + mg::outer( kb[i], edges[i] )) / scalarFactor;
//    o_plusGKB[i] = (2.0 * edgeMatrix[i - 1] - mg::outer( kb[i], edges[i - 1] )) / scalarFactor;
//    o_eqGKB[i] = -(o_minusGKB[i] + o_plusGKB[i]);

//}

