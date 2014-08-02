#include "ElasticRod.h"
#include "Utils.h"


ElasticRod::ElasticRod()
{
    m_bendStiffnes = 0.6;
    m_twistStiffness = 1.0;
    m_nIterations = 3;
}

ElasticRod::~ElasticRod()
{ }


void ElasticRod::init(const std::vector<mg::Vec4D>& restpos,
          const mg::Vec4D u0,
          const std::vector<mg::Vec4D>& pos,
          const std::vector<mg::Vec4D>& vel,
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


    m_u0 = u0;
    m_ppos = pos;
    m_pvel = vel;
    m_pmass = mass;
    m_twistAngle = twistAngle;
    m_pIsFixed = isFixed;

    m_totalTwist = 8.;

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
}


void ElasticRod::computeEdges(const std::vector<mg::Vec4D>& vertices,
                              std::vector<mg::Vec4D>& o_edges) const
{
    for (unsigned i = 0; i < vertices.size() - 1; ++i)
    {
        o_edges[i] = vertices[i + 1] - vertices[i];
    }
}

void ElasticRod::computeLengths(const std::vector<mg::Vec4D>& edges,
                                std::vector<mg::Real> &o_edgeL,
                                std::vector<mg::Real> &o_regionL,
                                mg::Real& o_totalL) const
{
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        o_edgeL[i] = edges[i].length();
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

void ElasticRod::computeKB(const std::vector<mg::Vec4D>& edges,
                           std::vector<mg::Vec4D>& o_kb) const
{
    o_kb[0].set(0, 0, 0, 0);

    double magnitude;
    double sinPhi, cosPhi;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
///     NOTE: as metioned in the paper the following formula produces |kb| = 2 * tan(phi / 2),
///     where phi is the angle of rotation defined by the 2 edges
        mg::Vec3D e1, e2;
        e1.set(edges[i - 1][0], edges[i - 1][1], edges[i - 1][2]);
        e2.set(edges[i][0], edges[i][1], edges[i][2]);
        e1 = 2 * cml::cross(e1, e2) /
                (m_restEdgeL[i - 1] * m_restEdgeL[i] + dot(edges[i - 1], edges[i]));

        o_kb[i].set(e1[0], e1[1], e1[2], 0);
//        o_kb[i] = 2 * cml::cross(mg::Vec3D(edges[i - 1]), mg::Vec3D(edges[i])) /
//                (m_restEdgeL[i - 1] * m_restEdgeL[i] + dot(edges[i - 1], edges[i]));

//      need to assert that |kb| = 2 * tan(phi / 2)
        magnitude = length_squared(o_kb[i]);
        extractSinAndCos(magnitude, sinPhi, cosPhi);
    }
}

void ElasticRod::extractSinAndCos(const double &magnitude,
                                  double &o_sinPhi, double &o_cosPhi) const
{
    o_cosPhi = sqrt(4.0 / (4.0 + magnitude));
    o_sinPhi = sqrt(magnitude / (4.0 + magnitude));

    assert( fabs((o_sinPhi * o_sinPhi + o_cosPhi * o_cosPhi) - 1.0) < mg::ERR );
    assert( fabs(magnitude - 4.0 * o_sinPhi * o_sinPhi / (o_cosPhi * o_cosPhi)) < mg::ERR );
}

void ElasticRod::computeBishopFrame(const mg::Vec4D& u0,
                        const std::vector<mg::Vec4D>& edges,
                        std::vector<mg::Vec4D>& o_kb,
                        std::vector<mg::Vec4D>& o_u,
                        std::vector<mg::Vec4D>& o_v) const
{
//    check if Bishop frame is correctly defined
    assert( fabs(length_squared(u0) - 1) < mg::ERR );
    assert( fabs(dot(u0, edges[0])) < mg::ERR );

    computeKB(edges, o_kb);

    o_u[0] = u0;
    mg::Vec3D e, u;
    e.set(edges[0][0], edges[0][1], edges[0][2]);
    u.set(o_u[0][0], o_u[0][1], o_u[0][2]);
    u = cml::cross(e, u);

    o_v[0].set(u[0], u[1], u[2], 0);
    o_v[0].normalize();

    double magnitude;
    double sinPhi, cosPhi;
//    compute Bishop frame for current configuration given u0 and edges
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        magnitude = length_squared(o_kb[i]);
        if (magnitude < mg::ERR)
        {
            o_u[i] = o_u[i - 1];
            o_v[i] = o_v[i - 1];
            continue;
        }
//        here sinPhi and cosPhi are derived from the length of kb |kb| = 2 * tan( phi/2 )
        extractSinAndCos(magnitude, sinPhi, cosPhi);
//        rotate frame u axis around kb
        mg::Quaternion q(cosPhi, sinPhi * normalize(o_kb[i]));
        mg::Quaternion p(0, o_u[i - 1][0], o_u[i - 1][1], o_u[i - 1][2]);
        p = q * p * conjugate(q);

        o_u[i].set(p[1], p[2], p[3], 0);

        e.set(edges[i][0], edges[i][1], edges[i][2]);
        u.set(o_u[i][0], o_u[i][1], o_u[i][2]);
        u = cml::cross(e, u);

        o_v[i].set(u[0], u[1], u[2], 0);
        o_v[i].normalize();

        assert( fabs(p[0]) <  mg::ERR );
        assert( fabs(dot(o_u[i], edges[i])) < mg::ERR );
        assert( fabs(length_squared(o_u[i]) - 1) < mg::ERR );
    }
}

void ElasticRod::computeMaterialFrame(const mg::Vec4D& u0,
                          const std::vector<mg::Vec4D>& edges,
                          const std::vector<mg::Real>& twistAngle,
                          std::vector<mg::Vec4D> &o_kb,
                          std::vector<mg::Vec4D>& o_m1,
                          std::vector<mg::Vec4D>& o_m2) const
{
    assert( edges.size() == (m_ppos.size() - 1) );
    assert( twistAngle.size() == edges.size() );

    computeBishopFrame(u0, edges, o_kb, o_m1, o_m2);

    mg::Real sinPhi, cosPhi;
    mg::Vec4D m1, m2;
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        sinPhi = sin(twistAngle[i]);
        cosPhi = cos(twistAngle[i]);
        m1 = cosPhi * o_m1[i] + sinPhi * o_m2[i];
        m2 = -sinPhi * o_m1[i] + cosPhi * o_m2[i];

        o_m1[i] = m1;
        o_m2[i] = m2;

        assert( fabs( dot(o_m1[i], o_m2[i])) < mg::ERR );
        assert( fabs( dot(o_m1[i], edges[i])) < mg::ERR );
        assert( fabs( dot(o_m2[i], edges[i])) < mg::ERR );
    }
}

void ElasticRod::computeMaterialCurvature(const std::vector<mg::Vec4D>& kb,
                                        const std::vector<mg::Vec4D>& m1,
                                        const std::vector<mg::Vec4D>& m2,
                                        std::vector<mg::Vec2D>& o_Wprev,
                                        std::vector<mg::Vec2D>& o_Wnext) const
{
    o_Wprev[0].set(0, 0);
    o_Wnext[0].set(0, 0);
    for (unsigned i = 1; i < kb.size(); ++i)
    {
        o_Wprev[i].set( dot( kb[i], m2[i - 1] ), -dot( kb[i], m1[i - 1] ));
        o_Wnext[i].set( dot( kb[i], m2[i] )    , -dot( kb[i], m1[i] )    );
    }
}

/* Use Runge-Kutta integration. */
//void HairStrand::update(mg::Real dt)
//{
//    std::vector<mg::Vec4D> k1, k2, k3, k4;
//    std::vector<mg::Vec4D> vertexScratch;

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

/* Use Verlet integration. */
void ElasticRod::update(mg::Real dt)
{
//    integrate centerline
    std::vector<mg::Vec4D> forces(m_ppos.size(), mg::Vec4D(0,0,0,0));
    computeForces(m_ppos, forces);

    std::vector<mg::Vec4D> prevPos(m_ppos.size());
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        prevPos[i] = m_ppos[i];
        m_pvel[i] += forces[i] * dt;
        m_ppos[i] += m_pvel[i] * dt;
    }

//    solve constraints
//    need FIX: check which ponts on the rod are with fixed position
    for (unsigned k = 0; k < m_nIterations; ++k)
    {
        for (unsigned i = 1; i < m_ppos.size(); ++i)
        {
            mg::Vec4D e = m_ppos[i - 1] - m_ppos[i];
            mg::Real l = e.length();
            if (l > mg::ERR)
            {
                e.normalize();
            }
            l = (l - m_restEdgeL[i - 1]);
            if (i == 1)
            {
                m_ppos[i] += l * e;
//            } else if (i == m_vertices.size() - 1)
//            {
//                m_vertices[i - 1] += -1.0 * l * e;
            }
            else {
                m_ppos[i] += 0.5 * l * e;
                m_ppos[i - 1] += -0.5 * l * e;
            }
        }
    }

//    update velocities:
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        m_pvel[i] = (m_ppos[i] - prevPos[i]) / dt;
    }
}

void ElasticRod::integrate(mg::Real dt)
{
//TODO: implement
}

void ElasticRod::solveConstraints(mg::Real dt)
{
//TODO: implement
}

void ElasticRod::computeForces(const std::vector<mg::Vec4D>& vertices, std::vector<mg::Vec4D>& o_forces)
{
//    computeExternalForces(vertices, o_forces);
//    computeElasticForces(vertices, o_forces);

    computeEdges(vertices, m_edges);
    computeKB(m_edges, m_kb);

    std::vector<mg::Vec4D> bendForces(vertices.size(), mg::Vec4D(0,0,0,0));
    std::vector<mg::Vec4D> twistForces(vertices.size(), mg::Vec4D(0,0,0,0));

    computeBendForces(vertices, m_edges, m_kb, bendForces);
    computeTwistForces(vertices, m_edges, m_kb, twistForces);

    for (unsigned i = 0; i < o_forces.size(); ++i)
    {
        if (!m_pIsFixed[i])
        {
            o_forces[i] += (bendForces[i] + twistForces[i] + mg::GRAVITY);//bendForces[i] + twistForces[i] + stretchScale * stretchForces[i]
        }
    }
}


void ElasticRod::computeExternalForces(const std::vector<mg::Vec4D>& vertices,
                                       std::vector<mg::Vec4D>& o_forces) const
{
    for (unsigned i = 1; i < o_forces.size(); ++i)
    {
        o_forces[i] += mg::GRAVITY;
    }
}


void ElasticRod::computeElasticForces(const std::vector<mg::Vec4D>& vertices,
                          std::vector<mg::Vec4D>& o_forces)
{
//    TODO: need to implement minimization of energy for twistAngles here

    computeEdges(m_ppos, m_edges);
    computeMaterialFrame(m_u0, m_edges, m_twistAngle, m_kb, m_m1, m_m2);
    computeMaterialCurvature(m_kb, m_m1, m_m2, m_Wprev,  m_Wnext);

    for (unsigned i = 1; i < o_forces.size(); ++i)
    {
        o_forces[i] += mg::GRAVITY;
    }
}


/**
 * Compute the forces acting on the vertices as a result of the
 * bending of the strand.
 **/
void ElasticRod::computeBendForces(const std::vector<mg::Vec4D>& vertices,
                                    const std::vector<mg::Vec4D>& edges,
                                    const std::vector<mg::Vec4D>& kb,
                                    std::vector<mg::Vec4D>& o_bendForces) const
{
    std::vector<mg::Matrix4D> edgeMatrices(edges.size());
    computeEdgeMatrices(edges, edgeMatrices);

    std::vector<mg::Vec4D> jminus_term(vertices.size(), mg::Vec4D(0,0,0,0));
    std::vector<mg::Vec4D> jplus_term(vertices.size(), mg::Vec4D(0,0,0,0));
    std::vector<mg::Vec4D> jequal_term(vertices.size(), mg::Vec4D(0,0,0,0));

    mg::Real scalarFactor;
    mg::Matrix4D minus_kbGradient, plus_kbGradient;
    for (unsigned int i = 1; i < o_bendForces.size() - 1; ++i)
    {
        scalarFactor = 1.0 / (m_restEdgeL[i - 1] * m_restEdgeL[i] + dot( edges[i - 1], edges[i]) );

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

void ElasticRod::computeTwistForces(const std::vector<mg::Vec4D>& vertices,
                                    const std::vector<mg::Vec4D>& edges,
                                    const std::vector<mg::Vec4D>& kb,
                                    std::vector<mg::Vec4D>& o_twistForces) const
{
    std::vector<mg::Vec4D> jminus_term(vertices.size(), mg::Vec4D(0,0,0,0));
    std::vector<mg::Vec4D> jplus_term(vertices.size(), mg::Vec4D(0,0,0,0));
    std::vector<mg::Vec4D> jequal_term(vertices.size(), mg::Vec4D(0,0,0,0));

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


/**
 * Computes skew-symmetric matrix 4x4 transpose( [e] ), such that [e] * x = cross( e, x )
 **/
void ElasticRod::computeEdgeMatrices(const std::vector<mg::Vec4D>& edges, std::vector<mg::Matrix4D>& o_edgeMatrices) const
{
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        mg::Matrix4D& em = o_edgeMatrices[i];
//        set 1st row
        em(0, 0) = 0;
        em(0, 1) = edges[i][2];
        em(0, 2) = -edges[i][1];
        em(0, 3) = 0;
//        set 2nd row
        em(1, 0) = -edges[i][2];
        em(1, 1) = 0;
        em(1, 2) = edges[i][0];
        em(1, 3) = 0;
//        set 3rd row
        em(2, 0) = edges[i][1];
        em(2, 1) = -edges[i][0];
        em(2, 2) = 0;
        em(2, 3) = 0;
//        set 4th column (translation)
        em(3, 0) = 0;
        em(3, 1) = 0;
        em(3, 2) = 0;
        em(3, 3) = 1;
    }
}
