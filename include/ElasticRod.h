#ifndef ROD_H
#define ROD_H

#include <vector>
#include "Types.h"

class ElasticRod {

public:

    ElasticRod();
    ~ElasticRod();

/// calculations are carried out in the coordinate space of pos
/// must assure everything is in same coordinate space
    void init(const std::vector<mg::Vec3D>& restpos,
              const mg::Vec3D u0,
              const std::vector<mg::Vec3D>& pos,
              const std::vector<mg::Vec3D>& vel,
              const std::vector<mg::Real>& mass,
              const std::vector<mg::Real> &twistAngle,
              const std::vector<bool>& isFixed
              );

    void update(mg::Real dt);

///    bending stiffness B = [EI1 0  ]   =  [bendStiffnes  0           ]
///                          [0   EI2]      [0             bendStiffnes]
///    depends on the Young modulus E(this is material property measured and constant throughout simulation)
///    and the area moment of inertia I1 or I2 which in turn depend only on cross-section area
///    note that I1 = I2 when the cross-section is circular, thus we always assume isotropic response
///    otherwise there should be 2 values describing the stiffness in each material axis direction
///    FIX: need to calculate bendStiffnes by using the thickness and density
    mg::Real m_bendStiffnes;
    mg::Matrix2D m_B;

///    twisting stiffness m_twistStiffness = GJ
///    depends on the shear modulus G = E/(2(1 + v)) where E is Young modulus and v is the Poison ratio
///    (this are material properties measured and constant throughout simulation)
///    and the area moment of twist J which in turn depends only on cross-section area
///    FIX: need to calculate twistStiffness by using the thickness and density
    mg::Real m_twistStiffness;

    unsigned m_nIterations;

    mg::Real m_maxElasticForceThreshold;

//    ========= should be private ========

/// #vertices = n + 1
    std::vector<mg::Vec3D> m_ppos;
/// #vel = n + 1
    std::vector<mg::Vec3D> m_pvel;
/// #mass = n + 1
    std::vector<mg::Real> m_pmass;
/// #isFixed = n + 1
    std::vector<bool> m_pIsFixed;
/// unit length vector - defines the Bishop frame of the first edge
    mg::Vec3D m_u0;

private:
/// #twistAngle = #edges = n
    std::vector<mg::Real> m_twistAngle;

/// #restEdgeL = #edges = n edge length
    std::vector<mg::Real> m_restEdgeL;
/// #restRegionL = #edges = n defines integrated length of voronoi region
/// NOTE: there is no restRegionL[0] => restRegionL[0] = 0
    std::vector<mg::Real> m_restRegionL;
/// #restWprev = #edges - restWprev[i] defines rest material curvatures at e[i -1]
    std::vector<mg::Vec2D> m_restWprev;
/// #restWnext = #edges - restWnext[i] defines rest material curvatures at e[i]
    std::vector<mg::Vec2D> m_restWnext;


/// #edges = #vertices - 1 defines tangent axis between v[i + 1] and v[i]
    std::vector<mg::Vec3D> m_edges;
/// defined at vertices 1 ... n - 1 i.e. #kb = (#edges - 1)
/// #kb = #edges - defines rotation from e[i-1] to e[i],
/// NOTE: there is no kb[0] => kb[0] = (0, 0, 0, 0);
    std::vector<mg::Vec3D> m_kb;
/// #m1 = #edges - m1[i] defines material axis m1 at e[i]
    std::vector<mg::Vec3D> m_m1;
/// #m2 = #edges - m2[i] defines material axis m2 at e[i]
    std::vector<mg::Vec3D> m_m2;

/// #restWprev = #edges - restWprev[i] defines rest material curvatures at e[i -1]
    std::vector<mg::Vec2D> m_Wprev;
/// #restWnext = #edges - restWnext[i] defines rest material curvatures at e[i]
    std::vector<mg::Vec2D> m_Wnext;


    mg::Real m_totalTwist;
    mg::Real m_totalL;




/// Computes the edge vectors between each consequent pair of vertices.
/// For n + 1 vertices there are n edges
    void computeEdges(const std::vector<mg::Vec3D>& vertices,
                      std::vector<mg::Vec3D>& o_edges) const;

/// Computes rest edge lengths the li scalars:
/// l[i] = |e[i - 1]| + |e[i + 1]|
/// NOTE: there is no l_0
    void computeLengths(const std::vector<mg::Vec3D> &edges,
                        std::vector<mg::Real> &o_edgeL,
                        std::vector<mg::Real> &o_regionL,
                        mg::Real &o_totalL) const;

    void extractSinAndCos(const double &magnitude,
                          double &o_sinPhi, double &o_cosPhi) const;

/// Computes the curvature binormal kb 3-vectors which depend on m_restEdgeL state
/// NOTE: there is no kb_0
/// depend on m_restEdgeL state
    void computeKB(const std::vector<mg::Vec3D>& edges,
                   std::vector<mg::Vec3D>& o_kb) const;

/// Computes Bishop frame for every edge by parallel transporting around the curvature binormal kb
/// the method internally computes kb 3-vectors
    void computeBishopFrame(const mg::Vec3D& u0,
                            const std::vector<mg::Vec3D>& edges,
                            std::vector<mg::Vec3D>& o_kb,
                            std::vector<mg::Vec3D>& o_u,
                            std::vector<mg::Vec3D>& o_v) const;

/// Computes Material frame for every edge - defining the orientation of the edge
/// the method internally computes kb 3-vectors, Bishop frame and rotates Bishop frame with twist angle
    void computeMaterialFrame(const mg::Vec3D& u0,
                              const std::vector<mg::Vec3D>& edges,
                              const std::vector<mg::Real>& twistAngle,
                              std::vector<mg::Vec3D>& o_kb,
                              std::vector<mg::Vec3D>& o_m1,
                              std::vector<mg::Vec3D>& o_m2) const;

    void computeMaterialCurvature(const std::vector<mg::Vec3D> &kb,
                                  const std::vector<mg::Vec3D>& m1,
                                  const std::vector<mg::Vec3D>& m2,
                                  std::vector<mg::Vec2D>& o_Wprev,
                                  std::vector<mg::Vec2D>& o_Wnext) const;

    void computeGradientKB(const std::vector<mg::Vec3D> &kb,
                           const std::vector<mg::Vec3D> &edges,
                           std::vector<mg::Matrix3D>& o_minusGKB,
                           std::vector<mg::Matrix3D>& o_plusGKB,
                           std::vector<mg::Matrix3D>& o_eqGKB) const;

    void computeGradientHolonomy(const std::vector<mg::Vec3D> &kb,
                           std::vector<mg::Vec3D>& o_minusGH,
                           std::vector<mg::Vec3D>& o_plusGH,
                           std::vector<mg::Vec3D>& o_eqGH) const;

    void computeGradientCurvature(int i, int k, int j,
                                const std::vector<mg::Matrix3D>& minusGKB,
                                const std::vector<mg::Matrix3D>& plusGKB,
                                const std::vector<mg::Matrix3D>& eqGKB,
                                const std::vector<mg::Vec3D>& minusGH,
                                const std::vector<mg::Vec3D>& plusGH,
                                const std::vector<mg::Vec3D>& eqGH,
                                const std::vector<mg::Vec2D>& wj,
                                mg::Matrix23D &o_GW) const;

/// Computes skew-symmetric matrix 4x4 transpose( [e] ), such that [e] * x = cross( e, x )
    void computeEdgeMatrices(const std::vector<mg::Vec3D>& edges,
                             std::vector<mg::Matrix3D>& o_edgeMatrices) const;


    void integrate(mg::Real dt);

    void solveConstraints(mg::Real dt);

    void computeForces(const std::vector<mg::Vec3D>& vertices,
                       std::vector<mg::Vec3D>& o_forces);

    void computeExternalForces(const std::vector<mg::Vec3D>& vertices,
                                std::vector<mg::Vec3D>& o_forces) const;

    void computeElasticForces(const std::vector<mg::Vec3D>& vertices,
                              std::vector<mg::Vec3D>& o_forces);

    void parallelTransportFrame(const mg::Vec3D& e0, const mg::Vec3D& e1,
                                mg::Vec3D& io_u) const;
/// Compute the forces acting on the vertices as a result of the
/// bending of the strand.
    void computeBendForces(const std::vector<mg::Vec3D>& vertices,
                            const std::vector<mg::Vec3D>& edges,
                            const std::vector<mg::Vec3D>& kb,
                            std::vector<mg::Vec3D>& o_bendForces) const;

/// Compute the forces acting on the vertices as a result of the
/// twisting of the strand.
    void computeTwistForces(const std::vector<mg::Vec3D>& vertices,
                            const std::vector<mg::Vec3D>& edges,
                            const std::vector<mg::Vec3D>& kb,
                            std::vector<mg::Vec3D>& o_twistForces) const;



};

#endif // ROD_H
