#ifndef ROD_H
#define ROD_H

#include <vector>
#include "ngl/Vec4.h"

class HairStrand {

public:

    HairStrand();
    ~HairStrand();

    void init();
    void update(ngl::Real dt);

///    bending stiffness B = [ EI1 0  ]   =  [ materialAlpha  0             ]
///                          [ 0   EI2]      [ 0              material alpha]
///    depends on the Young modulus E(this is material property measured and constant throughout simulation)
///    and the area moment of inertia I which in turn depends only on cross-section
///    FIX: need to calculate alpha by using the thickness and density
    ngl::Real m_bendStiffnes;
///    twisting stiffness m_twistStiffness = GJ
//    depends on the shear modulus G = E/(2(1 + v)) where v is the Poison ratio (this is material property measured and constant throughout simulation)
//    and the area moment of twist J which in turn depends only on cross-section
//    FIX: need to calculate alpha by using the thickness and density
    ngl::Real m_twistStiffness;

    ngl::Real m_length;
    ngl::Real m_density;
    ngl::Real m_thickness;
    unsigned m_nParticles;
    unsigned m_nIterations;

// =============== should be private ==================
    std::vector<ngl::Vec4> m_ppos;
    std::vector<ngl::Vec4> m_pvel;
    std::vector<bool> m_pIsFixed;

    std::vector<ngl::Vec4> m_edgesInit;
    std::vector<ngl::Real> m_l;
    ngl::Real m_totalL;

    std::vector<ngl::Vec4> m_verticesInit;
    std::vector<ngl::Vec4> m_velocitiesInit;
    ngl::Real m_totalTwist;

private:

//    note: the following are pure functions and don't change state
    void computeEdges(const std::vector<ngl::Vec4>& vertices, std::vector<ngl::Vec4>& o_edges) const;
    void computeLengths(const std::vector<ngl::Vec4> &edges, std::vector<ngl::Real> &o_edgeLength, ngl::Real &o_totalLength) const;

//    note: depend on current state
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

    void computeKB(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Vec4>& o_kb) const;
    void computeEdgeMatrices(const std::vector<ngl::Vec4>& edges, std::vector<ngl::Mat4>& o_edgeMatrices) const;
    void computeMatrixMult(const ngl::Vec4& kb, const ngl::Vec4& e, ngl::Mat4& o_matrix) const;

};

#endif // ROD_H
