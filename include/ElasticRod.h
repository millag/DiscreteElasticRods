#ifndef ROD_H
#define ROD_H

#include <vector>
#include "ngl/Vec4.h"

class ElasticRod {

public:

    ElasticRod();
    ~ElasticRod();

/// calculations are carried out in the coordinate space of pos
/// so user must assure everything is specified in the same coordinate space
    void init(const std::vector<ngl::Vec4>& restpos,
              const ngl::Vec4 u0,
              const std::vector<ngl::Vec4>& pos,
              const std::vector<ngl::Vec4>& vel,
              const std::vector<ngl::Real>& mass,
              const std::vector<ngl::Real> &twistAngle,
              const std::vector<bool>& isFixed
              );

    void update(ngl::Real dt);

///    bending stiffness B = [ EI1 0  ]   =  [ materialAlpha  0             ]
///                          [ 0   EI2]      [ 0              material alpha]
///    depends on the Young modulus E(this is material property measured and constant throughout simulation)
///    and the area moment of inertia I1 or I2 which in turn depend only on cross-section
///    note that I1 = I2 when the cross-section is circular, thus we always assume circular cross-section here
///    otherwise there should be 2 values describing the stiffness in each direction
///    FIX: need to calculate alpha by using the thickness and density
    ngl::Real m_bendStiffnes;
///    twisting stiffness m_twistStiffness = GJ
//    depends on the shear modulus G = E/(2(1 + v)) where v is the Poison ratio (this is material property measured and constant throughout simulation)
//    and the area moment of twist J which in turn depends only on cross-section
//    FIX: need to calculate alpha by using the thickness and density
    ngl::Real m_twistStiffness;
    unsigned m_nIterations;

//    ========= should be private ========

/// #vertices = n + 1
    std::vector<ngl::Vec4> m_ppos;
/// #vel = n + 1
    std::vector<ngl::Vec4> m_pvel;
/// #mass = n + 1
    std::vector<ngl::Real> m_pmass;
/// #isFixed = n + 1
    std::vector<bool> m_pIsFixed;
/// should be unit length and perpendicular to the first edge
    ngl::Vec4 m_u0;
/// #twistAngle = #edges = n
    std::vector<ngl::Real> m_twistAngle;


private:
/// #edges = #vertices - 1
    std::vector<ngl::Vec4> m_edges;
/// defined at vertices 1 ... n - 1 i.e. #kb = (#edges - 1)
/// #kb = #edges - defines rotation from e[i-1] to e[i], kb[0] = (0, 0, 0, 0);
    std::vector<ngl::Vec4> m_kb;
/// #m1 = #edges - m1[i] defines material axis m1 at e[i]
    std::vector<ngl::Vec4> m_m1;
/// #m2 = #edges - m2[i] defines material axis m2 at e[i]
    std::vector<ngl::Vec4> m_m2;

    std::vector<ngl::Real> m_restEdgeL;
    std::vector<ngl::Real> m_restRegionL;

    ngl::Real m_totalTwist;
    ngl::Real m_totalL;



//    void computeMaterialCurvature(const std::vector<ngl::Vec4>& kb, std::vector<ngl::Vec2>& o_w) const;


/// Computes the edge vectors between each consequent pair of vertices.
/// For n + 1 vertices there are n edges
    void computeEdges(const std::vector<ngl::Vec4>& vertices,
                      std::vector<ngl::Vec4>& o_edges) const;
/// Computes rest edge lengths the li scalars:
/// l[i] = |e[i - 1]| + |e[i + 1]|
/// NOTE: there is no l_0
    void computeLengths(const std::vector<ngl::Vec4> &edges,
                        std::vector<ngl::Real> &o_edgeL,
                        std::vector<ngl::Real> &o_regionL, ngl::Real &o_totalL) const;

/// Computes the curvature binormal kb 3-vectors
/// NOTE: there is no kb_0
/// depend on m_restEdgeL state
    void computeKB(const std::vector<ngl::Vec4>& edges,
                   std::vector<ngl::Vec4>& o_kb) const;

    void computeMaterialFrame(const ngl::Vec4& u0,
                              const std::vector<ngl::Vec4>& edges,
                              const std::vector<ngl::Real>& twistAngle,
                              std::vector<ngl::Vec4>& o_kb,
                              std::vector<ngl::Vec4>& o_m1,
                              std::vector<ngl::Vec4>& o_m2) const;

    void computeForces(const std::vector<ngl::Vec4>& vertices,
                       std::vector<ngl::Vec4>& o_forces) const;
    void computeBendForces(const std::vector<ngl::Vec4>& vertices,
                            const std::vector<ngl::Vec4>& edges,
                            const std::vector<ngl::Vec4>& kb,
                            std::vector<ngl::Vec4>& o_bendForces) const;
    void computeTwistForces(const std::vector<ngl::Vec4>& vertices,
                            const std::vector<ngl::Vec4>& edges,
                            const std::vector<ngl::Vec4>& kb,
                            std::vector<ngl::Vec4>& o_twistForces) const;
    void computeStretchForces(const std::vector<ngl::Vec4>& vertices,
                                const std::vector<ngl::Vec4>& edges,
                                std::vector<ngl::Vec4>& o_stretchForces) const;


    void computeEdgeMatrices(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Mat4>& o_edgeMatrices) const;
    void computeMatrixMult(const ngl::Vec4& kb, const ngl::Vec4& e, ngl::Mat4& o_matrix) const;

};

#endif // ROD_H
