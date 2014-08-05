#include "ElasticRod.h"

#include <cmath>

#include "Utils.h"


ElasticRod::ElasticRod()
{
    m_bendStiffnes = 0.7;
    m_twistStiffness = 1.0;
    m_nIterations = 4;

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
    m_Wprev.resize(m_edges.size());
    m_Wnext.resize(m_edges.size());

//    compute edges & lengths for rest shape
    computeEdges(restpos, m_edges);
    computeLengths(m_edges, m_restEdgeL, m_restRegionL, m_totalL);
//    compute kb and material frame for rest shape
    computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);
//    precompute material curvature for rest shape
    computeMaterialCurvature(m_kb, m_m1, m_m2, m_restWprev,  m_restWnext);

    computeEdges(m_ppos, m_edges);
    computeMaterialFrame(m_u0, m_edges, m_twistAngle, m_kb, m_m1, m_m2);
    computeMaterialCurvature(m_kb, m_m1, m_m2, m_Wprev,  m_Wnext);
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
                                std::vector<mg::Real> &o_regionL,
                                mg::Real& o_totalL) const
{
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        o_edgeL[i] = edges[i].length();
        assert( std::isfinite(o_edgeL[i]) && o_edgeL[i] > mg::ERR );
    }

    o_regionL[0] = 0.0;
    o_totalL = 0.0;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_regionL[i] = o_edgeL[i - 1] + o_edgeL[i];
        o_totalL += o_regionL[i];
    }
    o_totalL *= 0.5;
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


//      need to assert that |kb| = 2 * tan(phi / 2)
        magnitude = mg::length_squared(o_kb[i]);
        assert( std::isfinite(magnitude) );
        extractSinAndCos(magnitude, sinPhi, cosPhi);
    }
}

void ElasticRod::extractSinAndCos(const double &magnitude,
                                  double &o_sinPhi, double &o_cosPhi) const
{
    assert( std::isfinite(magnitude) && magnitude >= 0.0 );

    o_cosPhi = mg::sqrt_safe(4.0 / (4.0 + magnitude));
    o_sinPhi = mg::sqrt_safe(magnitude / (4.0 + magnitude));

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
    assert( fabs(mg::length_squared( u0 ) - 1) < mg::THRESHODL );
//    assert( fabs(mg::dot( u0, mg::normalize(edges[0]) )) < mg::ERR );

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
        assert( std::isfinite(magnitude) );

        mg::Real dotp = mg::dot(mg::normalize(m_edges[i - 1]), mg::normalize(m_edges[i]));
        if (fabs(dotp - 1) < mg::ERR )
        {
            o_u[i] = o_u[i - 1];
            o_v[i] = o_v[i - 1];
            continue;
        }
        if (fabs(dotp + 1) < mg::ERR)
        {
            o_u[i] = -o_u[i - 1];
            o_v[i] = -o_v[i - 1];
            continue;
        }

//        here sinPhi and cosPhi are derived from the length of kb |kb| = 2 * tan( phi/2 )
        extractSinAndCos(magnitude, sinPhi, cosPhi);
        assert( mg::length_squared(o_kb[i]) > mg::ERR );

//        rotate frame u axis around kb
        mg::Quaternion q(cosPhi, sinPhi * mg::normalize(o_kb[i]));
        mg::Quaternion p(0, o_u[i - 1]);
        p = q * p * mg::conjugate(q);

        o_u[i].set(p[1], p[2], p[3]);
        o_v[i] = mg::cross(edges[i], o_u[i]);
        o_v[i].normalize();

        assert( fabs(p[0]) <  mg::ERR );
        assert( fabs(mg::length_squared(o_u[i]) - 1) < mg::THRESHODL );
//        std::cout << "LENGTH u[" << i << "] " << mg::length_squared(o_u[i]) << std::endl;
//        std::cout << "ANGLE u edge[" << i << "] " << mg::deg( std::acos( mg::dot( o_u[i], mg::normalize(edges[i]) ) ) ) << std::endl;
//        std::cout << "ANGLE v edge[" << i << "] " << mg::deg( std::acos( mg::dot( o_v[i], mg::normalize(edges[i]) ) ) ) << std::endl;
//        assert( fabs(mg::dot( o_u[i], mg::normalize(edges[i]) )) < mg::ERR );
//        assert( fabs(mg::dot( o_v[i], mg::normalize(edges[i]) )) < mg::ERR );
    }
}

void ElasticRod::computeMaterialFrame(const mg::Vec3D& u0,
                          const std::vector<mg::Vec3D>& edges,
                          const std::vector<mg::Real>& twistAngle,
                          std::vector<mg::Vec3D> &o_kb,
                          std::vector<mg::Vec3D>& o_m1,
                          std::vector<mg::Vec3D>& o_m2) const
{
    assert( edges.size() == (m_ppos.size() - 1) );
    assert( twistAngle.size() == edges.size() );

    computeBishopFrame(u0, edges, o_kb, o_m1, o_m2);

    mg::Real sinPhi, cosPhi;
    mg::Vec3D m1, m2;
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        cosPhi = std::cos(twistAngle[i]);
        sinPhi = std::sqrt(1 - cosPhi * cosPhi);

        m1 = cosPhi * o_m1[i] + sinPhi * o_m2[i];
        m2 = -sinPhi * o_m1[i] + cosPhi * o_m2[i];

        o_m1[i] = m1;
        o_m2[i] = m2;

        assert( fabs( mg::length_squared( m1 ) - 1 ) < mg::THRESHODL );
        assert( fabs( mg::length_squared( m2 ) - 1 ) < mg::THRESHODL );
        assert( fabs( mg::dot( m1, m2 )) < mg::ERR );
//        std::cout << "ANGLE m1 edge[" << i << "] " << mg::deg( std::acos( mg::dot( m1, mg::normalize(edges[i]) ) ) ) << std::endl;
//        std::cout << "ANGLE m2 edge[" << i << "] " << mg::deg( std::acos( mg::dot( m2, mg::normalize(edges[i]) ) ) ) << std::endl;
//        assert( fabs( mg::dot( m1, mg::normalize(edges[i]) )) < mg::THRESHODL );
//        assert( fabs( mg::dot( m2, mg::normalize(edges[i]) )) < mg::THRESHODL );
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
        o_Wprev[i].set( mg::dot( kb[i], m2[i - 1] ), -mg::dot( kb[i], m1[i - 1] ));
        o_Wnext[i].set( mg::dot( kb[i], m2[i] )    , -mg::dot( kb[i], m1[i] )    );
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
        m_pvel[i] += forces[i] * dt;
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
                l1 = m_pmass[i] / (m_pmass[i] + m_pmass[i + 1]) * l;
                l2 = -m_pmass[i + 1] / (m_pmass[i] + m_pmass[i + 1]) * l;
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

//    parallel transport first frame in time
    mg::Vec3D t0 = mg::normalize(m_edges[0]);
    computeEdges(m_ppos, m_edges);
    mg::Vec3D t1 = mg::normalize(m_edges[0]);
    parallelTransportFrame(t0, t1, m_u0);

//    std::cout << "ANGLE u0 edge[0] " << mg::deg( std::acos( mg::dot(m_u0, t1) ) ) << std::endl;

    computeMaterialFrame(m_u0, m_edges, m_twistAngle, m_kb, m_m1, m_m2);
    computeMaterialCurvature(m_kb, m_m1, m_m2, m_Wprev,  m_Wnext);
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
        o_forces[i] += mg::GRAVITY;
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

    mg::Matrix23D gw;
    for (unsigned i = 0; i < vertices.size(); ++i)
    {
        if (m_pIsFixed[i])
        {
            continue;
        }

        for (unsigned k = 1; k < m_edges.size(); ++k)
        {

            computeGradientCurvature(i, k, k - 1,
                                     minusGKB, plusGKB, eqGKB,
                                     minusGH, plusGH, eqGH,
                                     m_Wprev,
                                     gw);
            o_forces[i] -= mg::transpose(gw) * (m_B * (m_Wprev[k] - m_restWprev[k]) / m_restRegionL[k]);

            assert( std::isfinite(mg::length_squared(o_forces[i])) );
            computeGradientCurvature(i, k, k,
                                     minusGKB, plusGKB, eqGKB,
                                     minusGH, plusGH, eqGH,
                                     m_Wnext,
                                     gw);
            o_forces[i] -= mg::transpose(gw) * (m_B * (m_Wnext[k] - m_restWnext[k]) / m_restRegionL[k]);

            assert( std::isfinite(mg::length_squared(o_forces[i])) );
        }

//        add dE/d(theta_n) FIX: ugly needs refactoring
        mg::Matrix2D J;
        mg::matrix_rotation_2D(J, mg::Constants::pi_over_2());
        unsigned n = m_edges.size() - 1;
        mg::Real mn = (m_twistAngle[n] - m_twistAngle[n - 1]);
        mg::Vec3D GH;
        GH.zero();
        if (i > 1)
        {
            GH += plusGH[i - 1];
        }
        if (i < m_edges.size() - 1)
        {
            GH += eqGH[i];
        }
        if (i < m_edges.size() - 2)
        {
            GH += minusGH[i + 1];
        }

        o_forces[i] += (mg::dot(m_Wnext[n], J * m_B * (m_Wnext[n] - m_restWnext[n])) +
                        2 * m_twistStiffness * mn / m_restRegionL[n]) * GH;

        mg::Real maxForce = 100;
        if (mg::length_squared(o_forces[i]) > maxForce * maxForce)
        {
            o_forces[i].normalize();
            o_forces[i] *= maxForce;
        }

        assert( std::isfinite(mg::length_squared(o_forces[i])) );
    }
}

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


void ElasticRod::computeGradientCurvature(int i, int k, int j,
                                        const std::vector<mg::Matrix3D>& minusGKB,
                                        const std::vector<mg::Matrix3D>& plusGKB,
                                        const std::vector<mg::Matrix3D>& eqGKB,
                                        const std::vector<mg::Vec3D>& minusGH,
                                        const std::vector<mg::Vec3D>& plusGH,
                                        const std::vector<mg::Vec3D>& eqGH,
                                        const std::vector<mg::Vec2D> &wj,
                                        mg::Matrix23D &o_GW) const
{
    assert(j == k || j == (k - 1));

//    need to make o_GW zero 3x2 matrix
    o_GW.zero();

    if ( (k < i - 1) || (k > i + 2) || (j > i + 1))
    {
        return;
    }
//    compute gradient KB(GKB) term
    if (k < i + 2)
    {
        mg::Vec3D m1 = m_m1[j];
        mg::Vec3D m2 = m_m2[j];
        o_GW(0,0) = m2[0];
        o_GW(0,1) = m2[1];
        o_GW(0,2) = m2[2];

        o_GW(1,0) = -m1[0];
        o_GW(1,1) = -m1[1];
        o_GW(1,2) = -m1[2];

        if (k == i - 1)
        {
            o_GW = o_GW * plusGKB[k];
        }
        if (k == i)
        {
            o_GW = o_GW * eqGKB[k];
        }
        if (k == i + 1)
        {
            o_GW = o_GW * minusGKB[k];
        }
    }
//    compute gradient Holonomy(GH) term
    mg::Matrix2D J;
    mg::matrix_rotation_2D(J, mg::Constants::pi_over_2());
    if (j == i - 1)
    {
        o_GW -= J * mg::outer(wj[k], plusGH[j]);
    }
    if (j == i)
    {
        o_GW -= J * mg::outer(wj[k], eqGH[j]);
    }
    if (j == i + 1)
    {
        o_GW -= J * mg::outer(wj[k], minusGH[j]);
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

void ElasticRod::parallelTransportFrame(const mg::Vec3D& t0, const mg::Vec3D& t1,
                                        mg::Vec3D& io_u) const
{
    mg::Vec3D axis = mg::cross(t0, t1);
    if (mg::length_squared(axis) > mg::ERR)
    {
        axis.normalize();
        mg::Quaternion q;
        mg::quaternion_rotation_axis_angle( q, axis, mg::acos_safe( mg::dot(t0, t1) ) );
        mg::Quaternion p(0, io_u);
        p = q * p * mg::conjugate(q);

        io_u.set(p[1], p[2], p[3]);
    }
    else
    {
        axis = mg::cross(t1, io_u);
        io_u = mg::cross(axis, t1);
    }
    if (mg::length_squared(io_u) < mg::ERR)
    {
        std::cout << "Something is really wrong!" << std::endl;
    }
    io_u.normalize();
}


////===================================== only for straight case ==============================================


//void ElasticRod::computeForces(const std::vector<mg::Vec3D>& vertices, std::vector<mg::Vec3D>& o_forces)
//{
//    computeEdges(vertices, m_edges);
//    computeKB(m_edges, m_kb);

//    std::vector<mg::Vec3D> bendForces(vertices.size(), mg::Vec3D(0,0,0));
//    std::vector<mg::Vec3D> twistForces(vertices.size(), mg::Vec3D(0,0,0));

//    computeBendForces(vertices, m_edges, m_kb, bendForces);
//    computeTwistForces(vertices, m_edges, m_kb, twistForces);

//    for (unsigned i = 0; i < o_forces.size(); ++i)
//    {
//        if (!m_pIsFixed[i])
//        {
//            o_forces[i] += (bendForces[i] + twistForces[i]) / m_pmass[i] + mg::GRAVITY;
//        }
//    }
//}

//void ElasticRod::computeEdgeMatrices(const std::vector<mg::Vec3D>& edges,
//                                     std::vector<mg::Matrix3D>& o_edgeMatrices) const
//{
//    for (unsigned i = 0; i < edges.size(); ++i)
//    {
//        mg::Matrix3D& em = o_edgeMatrices[i];
////        set 1st row
//        em(0, 0) = 0;
//        em(0, 1) = edges[i][2];
//        em(0, 2) = -edges[i][1];
////        set 2nd row
//        em(1, 0) = -edges[i][2];
//        em(1, 1) = 0;
//        em(1, 2) = edges[i][0];
////        set 3rd row
//        em(2, 0) = edges[i][1];
//        em(2, 1) = -edges[i][0];
//        em(2, 2) = 0;
//    }
//}

void ElasticRod::computeBendForces(const std::vector<mg::Vec3D>& vertices,
                                    const std::vector<mg::Vec3D>& edges,
                                    const std::vector<mg::Vec3D>& kb,
                                    std::vector<mg::Vec3D>& o_bendForces) const
{
    std::vector<mg::Matrix3D> edgeMatrices(edges.size());
    computeEdgeMatrices(edges, edgeMatrices);

    std::vector<mg::Vec3D> jminus_term(vertices.size(), mg::Vec3D(0,0,0));
    std::vector<mg::Vec3D> jplus_term(vertices.size(), mg::Vec3D(0,0,0));
    std::vector<mg::Vec3D> jequal_term(vertices.size(), mg::Vec3D(0,0,0));

    mg::Real scalarFactor;
    mg::Matrix3D minus_kbGradient, plus_kbGradient;
    for (unsigned int i = 1; i < o_bendForces.size() - 1; ++i)
    {
        scalarFactor = 1.0 / (m_restEdgeL[i - 1] * m_restEdgeL[i] + mg::dot( edges[i - 1], edges[i] ));

//        calculate  gradient of the curvature binormal kb for j = i - 1 and the corresponding force term
        minus_kbGradient = outer(edges[i], kb[i]);
        minus_kbGradient = (edgeMatrices[i] * 2.0 + minus_kbGradient) * scalarFactor;
        jminus_term[i] = minus_kbGradient * kb[i] * (-2 * m_bendStiffnes / m_restRegionL[i]);

//        calculate  gradient of the curvature binormal kb for j = i + 1 and the corresponding force term
        plus_kbGradient = outer(edges[i - 1], kb[i]);
        plus_kbGradient = (edgeMatrices[i - 1] * 2.0 + plus_kbGradient * (-1.0)) * scalarFactor;
        jplus_term[i] = plus_kbGradient * kb[i] * (-2 * m_bendStiffnes / m_restRegionL[i]);

//        calculate  gradient of the curvature binormal kb for j = i and the corresponding force term
        jequal_term[i] = ((minus_kbGradient + plus_kbGradient) * kb[i] * (-1.0)) * (-2 * m_bendStiffnes / m_restRegionL[i]);
    }

    for (unsigned i = 1; i < o_bendForces.size() - 1; ++i)
    {
        if (i > 1)
        {
            o_bendForces[i] += jplus_term[i - 1];
        }

        o_bendForces[i] += jequal_term[i];

        if (i < o_bendForces.size() - 2)
        {
            o_bendForces[i] += jminus_term[i + 1];
        }
    }

    o_bendForces[0] += jminus_term[1];
    o_bendForces[o_bendForces.size() - 1] += jplus_term[o_bendForces.size() - 2];
}

void ElasticRod::computeTwistForces(const std::vector<mg::Vec3D>& vertices,
                                    const std::vector<mg::Vec3D>& edges,
                                    const std::vector<mg::Vec3D>& kb,
                                    std::vector<mg::Vec3D>& o_twistForces) const
{
    std::vector<mg::Vec3D> jminus_term(vertices.size(), mg::Vec3D(0,0,0));
    std::vector<mg::Vec3D> jplus_term(vertices.size(), mg::Vec3D(0,0,0));
    std::vector<mg::Vec3D> jequal_term(vertices.size(), mg::Vec3D(0,0,0));

    double twistCoeff = m_twistStiffness * m_totalTwist / m_totalL;
    for (unsigned i = 1; i < o_twistForces.size() - 1; ++i)
    {
        jminus_term[i] = twistCoeff * kb[i] * 0.5 / m_restEdgeL[i - 1];
        jplus_term[i]  = twistCoeff * kb[i] * -0.5 / m_restEdgeL[i];
        jequal_term[i] = -(jminus_term[i] + jplus_term[i]);
    }

    for (unsigned i = 1; i < o_twistForces.size() - 1; ++i)
    {
        if (i > 1)
        {
            o_twistForces[i] += jplus_term[i - 1];
        }

        o_twistForces[i] += jequal_term[i];

        if (i < o_twistForces.size() - 2)
        {
            o_twistForces[i] += jminus_term[i + 1];
        }
    }

    o_twistForces[0] += jminus_term[1];
    o_twistForces[o_twistForces.size() - 1] += jplus_term[o_twistForces.size() - 2];
}


/* Use Runge-Kutta integration. */
//void HairStrand::update(mg::Real dt)
//{
//    std::vector<mg::Vec3D> k1, k2, k3, k4;
//    std::vector<mg::Vec3D> vertexScratch;

//    vertexScratch.resize(m_vertices.size());

//    computeForces(m_vertices, k1);
//    for (unsigned i = 0; i < m_vertices.size(); ++i)
//    {
//        vertexScratch[i] = m_vertices[i] + k1[i] * (0.5 * dt);
//    }

//    computeForces(vertexScratch, k2);
//    for (unsigned i = 0; i < m_vertices.size(); ++i)
//    {
//        vertexScratch[i] = m_vertices[i] + k2[i] * (0.5 * dt);
//    }

//    computeForces(vertexScratch, k3);
//    for (unsigned i = 0; i < m_vertices.size(); ++i)
//    {
//        vertexScratch[i] = m_vertices[i] + k3[i] * dt;
//    }

//    computeForces(vertexScratch, k4);
//    for (unsigned i = 0; i < m_vertices.size(); ++i)
//    {
//        //velocities[i] += forces[i] * timeStep;
//        //velocities[i] * timeStep;

//        m_vertices[i] += (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * dt / 6.0;
//    }
//}
