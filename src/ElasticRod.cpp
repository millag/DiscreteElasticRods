#include "ElasticRod.h"
#include "ngl/Mat4.h"
#include "ngl/Quaternion.h"
#include "Utils.h"



ElasticRod::ElasticRod()
{
    m_bendStiffnes = 0.6;
    m_twistStiffness = 1.0;


}

ElasticRod::~ElasticRod()
{ }


void ElasticRod::init(const std::vector<ngl::Vec4>& restpos,
          const ngl::Vec4 u0,
          const std::vector<ngl::Vec4>& pos,
          const std::vector<ngl::Vec4>& vel,
          const std::vector<ngl::Real>& mass,
          const std::vector<ngl::Real>& twistAngle,
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
//    init edges & lengths
    computeEdges(m_ppos, m_edges);

    computeLengths(m_edges, m_restEdgeL, m_restRegionL, m_totalL);
//    init kb and material frame
    computeMaterialFrame(m_u0, m_edges, m_twistAngle, m_kb, m_m1, m_m2);

    //TODO: precompute rest wi
    //compute material frame
}


void ElasticRod::computeEdges(const std::vector<ngl::Vec4>& vertices,
                              std::vector<ngl::Vec4>& o_edges) const
{
    o_edges.resize(vertices.size() - 1);
    for (unsigned i = 0; i < vertices.size() - 1; ++i)
    {
        o_edges[i] = vertices[i + 1] - vertices[i];
    }
}

void ElasticRod::computeLengths(const std::vector<ngl::Vec4>& edges,
                                std::vector<ngl::Real> &o_edgeL,
                                std::vector<ngl::Real> &o_regionL,
                                ngl::Real& o_totalL) const
{
    assert( edges.size() > 0 );

    o_edgeL.resize(edges.size());
    o_regionL.resize(edges.size());

    for (unsigned i = 0; i < o_edgeL.size(); ++i)
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

void ElasticRod::computeKB(const std::vector<ngl::Vec4>& edges,
                           std::vector<ngl::Vec4>& o_kb) const
{
    assert( edges.size() > 0 );
    assert( edges.size() == m_restEdgeL.size() );

    o_kb.resize(edges.size());
    o_kb[0].set(0, 0, 0, 0);

    ngl::Real magnitude, tan;
    ngl::Real sinPhi, cosPhi;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
///     NOTE: in the paper they claim the following formula produces |kb| = 2 * tan(phi / 2)
///     where phi is the angle between the 2 edges
        o_kb[i] = 2 * edges[i - 1].cross(edges[i]) /
                (m_restEdgeL[i - 1] * m_restEdgeL[i] + edges[i - 1].dot(edges[i]));

        assert( fabs(o_kb[i].dot(edges[i - 1])) < utils::ERR );
        assert( fabs(o_kb[i].dot(edges[i])) < utils::ERR );
//      need to assert that |kb| = 2 * tan(phi / 2)
        magnitude = o_kb[i].lengthSquared();
        if (magnitude  > 0.1)
            cosPhi = 0;
        cosPhi = sqrt(4.0 / (4.0 + magnitude));
        sinPhi = sqrt(magnitude / (4.0 + magnitude));
        tan = 4 * sinPhi * sinPhi / (cosPhi * cosPhi);
        assert( fabs((sinPhi * sinPhi + cosPhi * cosPhi) - 1) < utils::ERR );
        assert( fabs(magnitude - tan) < utils::ERR );
    }
}

void ElasticRod::computeMaterialFrame(const ngl::Vec4& u0,
                          const std::vector<ngl::Vec4>& edges,
                          const std::vector<ngl::Real>& twistAngle,
                          std::vector<ngl::Vec4> &o_kb,
                          std::vector<ngl::Vec4>& o_m1,
                          std::vector<ngl::Vec4>& o_m2) const
{
    assert( edges.size() == (m_ppos.size() - 1) );
    assert( twistAngle.size() == edges.size() );
//    check if Bishop frame is correctly defined
    assert( fabs(u0.lengthSquared() - 1) < utils::ERR );
    assert( fabs(u0.dot(edges[0])) < utils::ERR );

    o_m1.resize(edges.size());
    o_m2.resize(edges.size());

    computeKB(edges, o_kb);
    o_m1[0] = u0;
    o_m2[0] = edges[0].cross(o_m1[0]);
    o_m2[0].normalize();

    ngl::Real magnitude;
    ngl::Real sinPhi, cosPhi;
    ngl::Quaternion q, p;
//    compute Bishop frame for current configuration given u0 and edges
    for (unsigned i = 1; i < o_m1.size(); ++i)
    {
        magnitude = o_kb[i].lengthSquared();
        if (magnitude < utils::ERR)
        {
            o_m1[i] = o_m1[i - 1];
            o_m2[i] = o_m2[i - 1];
            continue;
        }
//        here sinPhi and cosPhi are derived from the length of kb |kb| = 2 * tan( phi/2 )
        cosPhi = sqrt(4.0 / (4.0 + magnitude));
        sinPhi = sqrt(magnitude / (4.0 + magnitude));
        assert( fabs(magnitude - 4 * sinPhi * sinPhi / (cosPhi * cosPhi)) < utils::ERR );

        magnitude = 1.0 / sqrt(magnitude);
        q.set(cosPhi, sinPhi * o_kb[i].m_x * magnitude, sinPhi * o_kb[i].m_y * magnitude, sinPhi * o_kb[i].m_z * magnitude);
        p.set(0, o_m1[i - 1].m_x, o_m1[i - 1].m_y, o_m1[i - 1].m_z);

        p = (q * p * q.conjugate());
        o_m1[i] = p.getVector();

        assert( fabs(p.getS()) <  utils::ERR );
        assert( fabs(o_m1[i].dot(edges[i])) < utils::ERR );
        assert( fabs(o_m1[i].lengthSquared() - 1) < utils::ERR );

        o_m2[i] = edges[i].cross(o_m1[i]);
        o_m2[i].normalize();
    }

    ngl::Vec4 m1, m2;
    for (unsigned i = 0; i < o_m1.size(); ++i)
    {
        sinPhi = sin(twistAngle[i]);
        cosPhi = cos(twistAngle[i]);

        m1 = cosPhi * o_m1[i] + sinPhi * o_m2[i];
        m2 = -sinPhi * o_m1[i] + cosPhi * o_m2[i];
        assert( fabs(m1.dot(m2)) < utils::ERR );
        assert( fabs(m1.dot(edges[i])) < utils::ERR );
        assert( fabs(m2.dot(edges[i])) < utils::ERR );

        o_m1[i] = m1;
        o_m2[i] = m2;
    }
}

/* Use Runge-Kutta integration. */
//void HairStrand::update(ngl::Real dt)
//{
//    std::vector<ngl::Vec4> k1, k2, k3, k4;
//    std::vector<ngl::Vec4> vertexScratch;

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
void ElasticRod::update(ngl::Real dt)
{
    std::vector<ngl::Vec4> forces;
    computeForces(m_ppos, forces);

    std::vector<ngl::Vec4> prevPos(m_ppos.size());
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        prevPos[i] = m_ppos[i];
        m_pvel[i] += forces[i] * dt;
        m_ppos[i] += m_pvel[i] * dt;
    }

    // solve constraints:
    // need FIX: check which ponts on the rod are with fixed position
    const unsigned nIter = 4;
    for (unsigned k = 0; k < nIter; ++k)
    {
        for (unsigned i = 1; i < m_ppos.size(); ++i)
        {
            ngl::Vec4 e = m_ppos[i - 1] - m_ppos[i];
            ngl::Real l = e.length();
            if (l > utils::ERR)
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

    // update velocities:
    for (unsigned i = 0; i < m_ppos.size(); ++i)
    {
        m_pvel[i] = (m_ppos[i] - prevPos[i]) / dt;
    }
}

void ElasticRod::computeForces(const std::vector<ngl::Vec4>& vertices, std::vector<ngl::Vec4>& o_forces) const
{
    o_forces.resize(vertices.size());

    std::vector<ngl::Vec4> edges;
    computeEdges(vertices, edges);
    std::vector<ngl::Vec4> kb;
    computeKB(edges, kb);

    std::vector<ngl::Vec4> bendForces(vertices.size());
    std::vector<ngl::Vec4> twistForces(vertices.size());

    computeBendForces(vertices, edges, kb, bendForces);
    computeTwistForces(vertices, edges, kb, twistForces);

    for (unsigned i = 0; i < o_forces.size(); ++i)
    {
        if (!m_pIsFixed[i])
        {
            o_forces[i] = (bendForces[i] + twistForces[i] + utils::G);//bendForces[i] + twistForces[i] + stretchScale * stretchForces[i]
        }
    }
}

/**
 * Compute the forces acting on the vertices as a result of the
 * bending of the strand.
 **/
void ElasticRod::computeBendForces(const std::vector<ngl::Vec4>& vertices,
                                    const std::vector<ngl::Vec4>& edges,
                                    const std::vector<ngl::Vec4>& kb,
                                    std::vector<ngl::Vec4>& o_bendForces) const
{
    o_bendForces.resize(vertices.size());

    std::vector<ngl::Mat4> edgeMatrices;
    computeEdgeMatrices(edges, edgeMatrices);

    std::vector<ngl::Vec4> jminus_term(vertices.size());
    std::vector<ngl::Vec4> jplus_term(vertices.size());
    std::vector<ngl::Vec4> jequal_term(vertices.size());

    ngl::Mat4 minus_kbGradient, plus_kbGradient;
    for (unsigned int i = 1; i < o_bendForces.size() - 1; ++i)
    {
        double scalarFactor = 1.0 / (m_restEdgeL[i - 1] * m_restEdgeL[i] + edges[i - 1].dot(edges[i]));

//        calculate  gradient of the curvature binormal kb for j = i - 1 and the corresponding force term
        computeMatrixMult(kb[i], edges[i], minus_kbGradient);
        minus_kbGradient = (edgeMatrices[i] * 2.0 + minus_kbGradient) * scalarFactor;
        jminus_term[i] = minus_kbGradient * kb[i] * (-2 * m_bendStiffnes / m_restRegionL[i]);

//        calculate  gradient of the curvature binormal kb for j = i + 1 and the corresponding force term
        computeMatrixMult(kb[i], edges[i - 1], plus_kbGradient);
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

    o_bendForces[0] = jminus_term[1];
    o_bendForces[o_bendForces.size() - 1] = jplus_term[o_bendForces.size() - 2];
}

void ElasticRod::computeTwistForces(const std::vector<ngl::Vec4>& vertices,
                                    const std::vector<ngl::Vec4>& edges,
                                    const std::vector<ngl::Vec4>& kb,
                                    std::vector<ngl::Vec4>& o_twistForces) const
{
    o_twistForces.resize(vertices.size());

    std::vector<ngl::Vec4> jminus_term(vertices.size());
    std::vector<ngl::Vec4> jplus_term(vertices.size());
    std::vector<ngl::Vec4> jequal_term(vertices.size());

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

    o_twistForces[0] = jminus_term[1];
    o_twistForces[o_twistForces.size() - 1] = jplus_term[o_twistForces.size() - 2];
}

/**
 * Computes the forces used to enforce the inextensibility constraint.
 **/
void ElasticRod::computeStretchForces(const std::vector<ngl::Vec4>& vertices,
                                        const std::vector<ngl::Vec4>& edges,
                                        std::vector<ngl::Vec4>& o_stretchForces) const
{
    o_stretchForces.resize(vertices.size());
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_stretchForces[i] = edges[i] * (4.0 / m_restEdgeL[i]) *
                            (edges[i].dot(edges[i]) - m_restEdgeL[i] * m_restEdgeL[i]);
        o_stretchForces[i] -= edges[i - 1] * (4.0 / m_restEdgeL[i - 1]) *
                            (edges[i - 1].dot(edges[i - 1]) - m_restEdgeL[i - 1] * m_restEdgeL[i - 1]);
    }

    unsigned idx = 0;
    o_stretchForces[idx] = edges[idx] * (4.0 / m_restEdgeL[idx]) *
                        (edges[idx].dot(edges[idx]) - m_restEdgeL[idx] * m_restEdgeL[idx]);

    idx = edges.size();

    assert(idx == (o_stretchForces.size() - 1));

    o_stretchForces[idx] = (-1.0) * edges[idx - 1] * (4.0 / m_restEdgeL[idx - 1]) *
                           (edges[idx - 1].dot(edges[idx - 1]) - m_restEdgeL[idx - 1] * m_restEdgeL[idx - 1]);

}


/**
 * Computes skew-symmetric matrix 4x4 [e], such that [e] * x = e.cross(x)
 **/
void ElasticRod::computeEdgeMatrices(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Mat4>& o_edgeMatrices) const
{
    assert(edges.size() == m_restEdgeL.size());

    o_edgeMatrices.resize(edges.size());

    for (unsigned i = 0; i < edges.size(); ++i)
    {
        ngl::Mat4& em = o_edgeMatrices[i];
//        set 1st column (xaxis)
        em.m_00 = 0;
        em.m_01 = edges[i].m_z;
        em.m_02 = -edges[i].m_y;
        em.m_03 = 0;
//        set 2nd column (yaxis)
        em.m_10 = -edges[i].m_z;
        em.m_11 = 0;
        em.m_12 = edges[i].m_x;
        em.m_13 = 0;
//        set 3rd column (zaxis)
        em.m_20 = edges[i].m_y;
        em.m_21 = -edges[i].m_x;
        em.m_22 = 0;
        em.m_23 = 0;
//        set 4th column (translation)
        em.m_30 = 0;
        em.m_31 = 0;
        em.m_32 = 0;
        em.m_33 = 1;
    }
}

void ElasticRod::computeMatrixMult(const ngl::Vec4& kb, const ngl::Vec4& e, ngl::Mat4& o_matrix) const
{
//        set 1st column (xaxis)
        o_matrix.m_00 = kb.m_x * e.m_x;
        o_matrix.m_01 = kb.m_y * e.m_x;
        o_matrix.m_02 = kb.m_z * e.m_x;
        o_matrix.m_03 = 0;
//        set 2nd column (yaxis)
        o_matrix.m_10 = kb.m_x * e.m_y;
        o_matrix.m_11 = kb.m_y * e.m_y;
        o_matrix.m_12 = kb.m_z * e.m_y;
        o_matrix.m_13 = 0;
//        set 3rd column (zaxis)
        o_matrix.m_20 = kb.m_x * e.m_z;
        o_matrix.m_21 = kb.m_y * e.m_z;
        o_matrix.m_22 = kb.m_z * e.m_z;
        o_matrix.m_23 = 0;
//        set 4th column (translation)
        o_matrix.m_30 = 0;
        o_matrix.m_31 = 0;
        o_matrix.m_32 = 0;
        o_matrix.m_33 = 1;
}
