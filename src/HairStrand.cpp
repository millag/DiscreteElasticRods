#include "HairStrand.h"
#include "ngl/Mat4.h"
#include "Utils.h"

HairStrand::HairStrand()
{
    m_bendStiffnes = 0.6;
    m_twistStiffness = 1.0;


}

HairStrand::~HairStrand()
{ }


void HairStrand::init()
{
//    init vertices, velocities, isFixed, totalTwist
    unsigned nVertices = 11;
    ngl::Vec4 start(0, 3, 0, 1);
    ngl::Vec4 end(3, 3, 0, 1);

    m_verticesInit.resize(nVertices);
    m_velocitiesInit.resize(nVertices);
    m_pIsFixed.resize(nVertices);

    ngl::Real t = 0.;
    for (unsigned i = 0; i < m_verticesInit.size(); ++i)
    {
        t = (ngl::Real)(i) / (nVertices - 1);
        m_verticesInit[i] = (1 - t) * start + t * end;
        m_velocitiesInit[i].set(0, 0, 0, 0);
        m_pIsFixed[i] = 0;
    }
    m_pIsFixed[0] = 1;
//    m_isFixedVertex[nVertices - 1] = 1;
    m_totalTwist = 8.;

//    reset
    m_ppos = m_verticesInit;
    m_pvel = m_velocitiesInit;

//    init edges & lengths
    computeEdges(m_verticesInit, m_edgesInit);
    computeLengths(m_edgesInit ,m_l, m_totalL);
}

/**
 * Computes the edge vectors between each consequent pair of vertices.
 * For n + 1 vertices there are n edges
 **/
void HairStrand::computeEdges(const std::vector<ngl::Vec4>& vertices, std::vector<ngl::Vec4>& o_edges) const
{
    o_edges.resize(vertices.size() - 1);
    for (unsigned i = 0; i < vertices.size() - 1; ++i)
    {
        o_edges[i] = vertices[i + 1] - vertices[i];
    }
}

/**
 * Computes the l scalars:
 *  l[i] = |e[i - 1]| + |e[i + 1]|
 *
 *  note: there is no l_0
 **/
void HairStrand::computeLengths(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Real>& o_edgeLength, ngl::Real& o_totalLength) const
{
    o_edgeLength.resize(edges.size());
    if (o_edgeLength.size() > 0)
    {
        o_edgeLength[0] = 0.0;
    }

    o_totalLength = 0.;
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_edgeLength[i] = edges[i - 1].length() + edges[i].length();
        o_totalLength += o_edgeLength[i];
    }

    o_totalLength *= 0.5;
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
void HairStrand::update(ngl::Real dt)
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
            l = (l - m_edgesInit[i - 1].length());
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

void HairStrand::computeForces(const std::vector<ngl::Vec4>& vertices, std::vector<ngl::Vec4>& o_forces) const
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
 * Computes the curvature binormal kb 3-vectors
 *
 * NOTE: there is no kb_0
 **/
void HairStrand::computeKB(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Vec4>& o_kb) const
{
    assert(edges.size() == m_edgesInit.size());

    o_kb.resize(edges.size());
    o_kb[0].set(0, 0, 0, 0);
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_kb[i] = 2 * edges[i - 1].cross(edges[i]) /
                (m_edgesInit[i - 1].length() * m_edgesInit[i].length() + edges[i - 1].dot(edges[i]));
    }
}

/**
 * Compute the forces acting on the vertices as a result of the
 * bending of the strand.
 **/
void HairStrand::computeBendForces(const std::vector<ngl::Vec4>& vertices,
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
        double scalarFactor = 1.0 / (m_edgesInit[i - 1].length() * m_edgesInit[i].length() + edges[i - 1].dot(edges[i]));

//        calculate  gradient of the curvature binormal kb for j = i - 1 and the corresponding force term
        computeMatrixMult(kb[i], edges[i], minus_kbGradient);
        minus_kbGradient = (edgeMatrices[i] * 2.0 + minus_kbGradient) * scalarFactor;
        jminus_term[i] = minus_kbGradient * kb[i] * (-2 * m_bendStiffnes / m_l[i]);

//        calculate  gradient of the curvature binormal kb for j = i + 1 and the corresponding force term
        computeMatrixMult(kb[i], edges[i - 1], plus_kbGradient);
        plus_kbGradient = (edgeMatrices[i - 1] * 2.0 + plus_kbGradient * (-1.0)) * scalarFactor;
        jplus_term[i] = plus_kbGradient * kb[i] * (-2 * m_bendStiffnes / m_l[i]);

//        calculate  gradient of the curvature binormal kb for j = i and the corresponding force term
        jequal_term[i] = ((minus_kbGradient + plus_kbGradient) * kb[i] * (-1.0)) * (-2 * m_bendStiffnes / m_l[i]);
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

void HairStrand::computeTwistForces(const std::vector<ngl::Vec4>& vertices,
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
        jminus_term[i] = twistCoeff * kb[i] * 0.5 / m_edgesInit[i - 1].length();
        jplus_term[i]  = twistCoeff * kb[i] * -0.5 / m_edgesInit[i].length();
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
void HairStrand::computeStretchForces(const std::vector<ngl::Vec4>& vertices,
                                        const std::vector<ngl::Vec4>& edges,
                                        std::vector<ngl::Vec4>& o_stretchForces) const
{
    o_stretchForces.resize(vertices.size());
    for (unsigned i = 1; i < edges.size(); ++i)
    {
        o_stretchForces[i] = edges[i] * (4.0 / m_edgesInit[i].length()) *
                            (edges[i].dot(edges[i]) - m_edgesInit[i].dot(m_edgesInit[i]));
        o_stretchForces[i] -= edges[i - 1] * (4.0 / m_edgesInit[i - 1].length()) *
                            (edges[i - 1].dot(edges[i - 1]) - m_edgesInit[i - 1].dot(m_edgesInit[i - 1]));
    }

    unsigned idx = 0;
    o_stretchForces[idx] = edges[idx] * (4.0 / m_edgesInit[idx].length()) *
                        (edges[idx].dot(edges[idx]) - m_edgesInit[idx].dot(m_edgesInit[idx]));

    idx = edges.size();

    assert(idx == (o_stretchForces.size() - 1));

    o_stretchForces[idx] = (-1.0) * edges[idx - 1] * (4.0 / m_edgesInit[idx - 1].length()) *
                           (edges[idx - 1].dot(edges[idx - 1]) - m_edgesInit[idx - 1].dot(m_edgesInit[idx - 1]));

}


/**
 * Computes skew-symmetric matrix 4x4 [e], such that [e] * x = e.cross(x)
 **/
void HairStrand::computeEdgeMatrices(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Mat4>& o_edgeMatrices) const
{
    assert(edges.size() == m_edgesInit.size());

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

void HairStrand::computeMatrixMult(const ngl::Vec4& kb, const ngl::Vec4& e, ngl::Mat4& o_matrix) const
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
