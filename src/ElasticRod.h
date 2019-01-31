#pragma once

#include "config.h"
#include <dlib/matrix.h>


typedef dlib::matrix<double,0,1> ColumnVector;

struct ElasticRodState
{
	friend class ElasticRod;
public:
	inline void clear()
	{
		m_ppos.clear();
		m_pvel.clear();
		m_u0.zero();
	}

private:
	std::vector<mg::Vec3D> m_ppos;
	std::vector<mg::Vec3D> m_pvel;
	mg::Vec3D m_u0;
};

struct ElasticRodParams
{
	enum class MinimizationStrategy
	{
		NONE = 0,
		NEWTON = 1,
		BFGS = 2,
		BFGS_NUMERIC = 3,
		SENTINEL,
	};

	ElasticRodParams( mg::Real bendStiffness = 1.0,
	          mg::Real twistStiffness = 1.0,
	          mg::Real maxElasticForce = 1000,
	          MinimizationStrategy strategy = MinimizationStrategy::BFGS,
	          double tolerance = 1e-6f,
	          unsigned maxIter = 100 ):
	    m_beta(twistStiffness)
	  , m_maxElasticForce(maxElasticForce)
	  , m_strategy(strategy)
	  , m_tolerance(tolerance)
	  , m_maxIter(maxIter)
	{
		setBendStiffness(bendStiffness);
	}

	inline void setBendStiffness(const mg::Real& bendStiffness)
	{
		m_B.identity();
		m_B *= bendStiffness;
	}

	inline void setTwistStiffness(const mg::Real& twistStiffness)
	{
		m_beta = twistStiffness;
	}

///    bending stiffness B = [EI1 0  ]   =  [bendStiffnes  0           ]
///                          [0   EI2]      [0             bendStiffnes]
///    depends on the Young modulus E(this is material property measured and constant throughout simulation)
///    and the area moment of inertia I1 or I2 which in turn depend only on cross-section area
///    note that I1 = I2 when the cross-section is circular, thus we always assume isotropic response
///    otherwise there should be 2 values describing the stiffness in each material axis direction
	mg::Matrix2D m_B;
///    twisting stiffness m_beta = GJ
///    depends on the shear modulus G = E/(2(1 + v)) where E is Young modulus and v is the Poison ratio
///    (this are material properties measured and constant throughout simulation)
///    and the area moment of twist J which in turn depends only on cross-section area
	mg::Real m_beta;
///    parameter that limits elastic force magnitude
///    need to limit the force otherwise when e[i] ~ -e[i - 1], ||kb|| goes to inf
///    => elastic force goes to inf and simulation blows up
	mg::Real m_maxElasticForce;
///     current twist angles drive the bending force for next frame
///     in order to decouple twist from bending energy minimization is required
///     elastic force computations relies on potential enery minimization of the rod
///     with respect to twist angle of the material frame
///     energy minimization is implemented using dlib numeric library
///     the following parameters control minimization method, tolerance and max iterartion used
///     consult dlib documentation for reference
	MinimizationStrategy m_strategy;
	double m_tolerance;
	unsigned m_maxIter;
};

class ElasticRod
{
public:
	ElasticRod();
	~ElasticRod();

	ElasticRod( ElasticRod&& other );
	ElasticRod& operator=( ElasticRod&& other );


/// calculations are carried out in the coordinate space of pos - you must assure everything is in same coordinate space
	void initialize(const ElasticRodParams& params,
	                const std::vector<mg::Vec3D>& restpos,
	                const mg::Vec3D& u0,
	                const std::vector<mg::Vec3D>& pos,
	                const std::vector<mg::Vec3D>& vel,
	                const std::vector<mg::Real>& mass,
	                const ColumnVector& theta,
	                const std::set<unsigned>& isClamped);

///    apply one PBD iteration for solving internal distance constraints
	void applyInternalConstraintsIteration();

///    calculate internal elastic forces for each position and accumulate the result in o_forces
///    first compute valid state based on current position and edge data
///    NOTE: the current state of the rod is not altered except when debug(DBUGG) is enabled
///    When DBUGG is on computed forces are stored for visualization purpose
	void accumulateInternalElasticForces(std::vector<mg::Vec3D>& o_forces);

	void getState(ElasticRodState& o_state) const;
	void setState(const ElasticRodState& state);

///    compute valid state based on current positions and edge data from previous frame
	void updateCurrentState();

public:
//    ========= must be private ========

/// #vertices = n + 1
	std::vector<mg::Vec3D> m_ppos;
/// #vel = n + 1
	std::vector<mg::Vec3D> m_pvel;
/// #mass = n + 1
	std::vector<mg::Real> m_pmass;
/// map containg indices of clamped positions - defines boundary conditions
	std::set<unsigned> m_isClamped;
/// unit length vector - defines the Bishop frame for the first edge for rest shape (restpos)
/// Bishop frame for the first edge for the current configuration(pos) is deduced by parallel transport
	mg::Vec3D m_u0;
/// #twistAngle = #edges = n
	ColumnVector m_theta;

/// defined at vertices 1 ... n - 1 i.e. #kb = (#edges - 1)
/// #kb = #edges - defines rotation from e[i-1] to e[i],
/// NOTE: there is no kb[0] => kb[0] = (0, 0, 0, 0);
	std::vector<mg::Vec3D> m_kb;
/// #m1 = #edges - m1[i] defines material axis m1(normal) at e[i]
/// NOTE: it is used also as temporary place holder for u axis of Bishop frame
	std::vector<mg::Vec3D> m_m1;
/// #m2 = #edges - m2[i] defines material axis m2(binormal) at e[i]
/// NOTE: it is used also as temporary place holder for v axis of Bishop frame
	std::vector<mg::Vec3D> m_m2;

/// #m_elasticForce = #positions - m_elasticForce[i] contains applied elastic force from last iteration
/// NOTE: this is used for debug purpose only - in general theres no need to store it
#ifdef DBUGG
	std::vector<mg::Vec3D> m_elasticForce;
#endif

private:
/// #restEdgeL = #edges = n edge length
	std::vector<mg::Real> m_restEdgeL;
/// #restRegionL = #edges = n defines integrated length for voronoi region of integration
/// NOTE: there is no restRegionL[0] => restRegionL[0] = 0
	std::vector<mg::Real> m_restRegionL;
/// #restWprev = #edges - restWprev[i] defines rest material curvature for kb[i] at e[i -1]
/// NOTE: there is no m_restWprev[0], because kb[0] = (0,0,0) => m_restWprev[0] = (0,0)
	std::vector<mg::Vec2D> m_restWprev;
/// #restWnext = #edges - restWnext[i] defines rest material curvature for kb[i] at e[i]
/// NOTE: there is no m_restWnext[0], because kb[0] = (0,0,0) => m_restWnext[0] = (0,0)
	std::vector<mg::Vec2D> m_restWnext;
/// #edges = #vertices - 1 defines tangent axis for the segment between v[i + 1] and v[i]
/// together with m_m1 and m_m2 they define material frame for each segment of the rod
	std::vector<mg::Vec3D> m_edges;

//    contains all params that drive the rod and can be changed dynamically (animated, whatever)
	const ElasticRodParams* m_params = nullptr;

	class MinimizationPImpl;
	std::unique_ptr<MinimizationPImpl> m_minimization;

private:
/// Computes the edge vectors between each consequent pair of vertices.
/// For n + 1 vertices there are n edges
	void computeEdges(const std::vector<mg::Vec3D>& vertices,
	                  std::vector<mg::Vec3D>& o_edges) const;

/// Computes rest edge lengths the li scalars:
/// l[i] = |e[i - 1]| + |e[i + 1]|
/// NOTE: there is no l[0]
	void computeLengths(const std::vector<mg::Vec3D>& vertices,
	                    std::vector<mg::Real>& o_edgeL,
	                    std::vector<mg::Real>& o_regionL)const;

	void computeW(const mg::Vec3D& kb,
	              const mg::Vec3D& m1,
	              const mg::Vec3D& m2,
	              mg::Vec2D& o_wij) const;

	void extractSinAndCos(const double& magnitude,
	                      double& o_sinPhi, double& o_cosPhi) const;

/// Computes the curvature binormal kb 3-vectors which depend on m_restEdgeL state
/// NOTE: there is no kb[0]
/// depend on m_restEdgeL state
	void computeKB(const std::vector<mg::Vec3D>& edges,
	               std::vector<mg::Vec3D>& o_kb) const;

/// Computes Bishop frame for every edge by parallel transporting u0 around the curvature binormal kb
/// the method internally computes kb 3-vectors
	void computeBishopFrame(const mg::Vec3D& u0,
	                        const std::vector<mg::Vec3D>& edges,
	                        std::vector<mg::Vec3D>& o_kb,
	                        std::vector<mg::Vec3D>& o_u,
	                        std::vector<mg::Vec3D>& o_v) const;

/// Computes Material frame for every edge - defining the orientation of the edge
/// rotates Bishop frame with twist angle theta
/// NOTE: the method expects io_m1 and io_m2 to contain u and v axis of the Bishop frame correspondingly
	void computeMaterialFrame(const ColumnVector& theta,
	                          std::vector<mg::Vec3D>& io_m1,
	                          std::vector<mg::Vec3D>& io_m2) const;

	void computeMaterialCurvature(const std::vector<mg::Vec3D>& kb,
	                              const std::vector<mg::Vec3D>& m1,
	                              const std::vector<mg::Vec3D>& m2,
	                              std::vector<mg::Vec2D>& o_Wprev,
	                              std::vector<mg::Vec2D>& o_Wnext) const;

	void parallelTransportFrame(const mg::Vec3D& e0,
	                            const mg::Vec3D& e1,
	                            mg::Vec3D& io_u) const;



	void computeGradientKB(const std::vector<mg::Vec3D>& kb,
	                       const std::vector<mg::Vec3D>& edges,
	                       std::vector<mg::Matrix3D>& o_minusGKB,
	                       std::vector<mg::Matrix3D>& o_plusGKB,
	                       std::vector<mg::Matrix3D>& o_eqGKB) const;

	void computeGradientHolonomyTerms(const std::vector<mg::Vec3D>& kb,
	                                   std::vector<mg::Vec3D>& o_minusGH,
	                                   std::vector<mg::Vec3D>& o_plusGH,
	                                   std::vector<mg::Vec3D>& o_eqGH) const;

	void computeGradientHolonomy(unsigned i , unsigned j,
	                             const std::vector<mg::Vec3D>& minusGH,
	                             const std::vector<mg::Vec3D>& plusGH,
	                             const std::vector<mg::Vec3D>& eqGH,
	                             mg::Vec3D& o_GH) const;

	void computeGradientCurvature(unsigned i, unsigned k, unsigned j,
	                            const std::vector<mg::Matrix3D>& minusGKB,
	                            const std::vector<mg::Matrix3D>& plusGKB,
	                            const std::vector<mg::Matrix3D>& eqGKB,
	                            const std::vector<mg::Vec3D>& minusGH,
	                            const std::vector<mg::Vec3D>& plusGH,
	                            const std::vector<mg::Vec3D>& eqGH,
	                            const mg::Vec2D& wkj,
	                            const mg::Matrix2D& J,
	                            mg::Matrix23D& o_GW) const;

	void computeEnergy(const std::vector<mg::Vec3D>& m1,
	                   const std::vector<mg::Vec3D>& m2,
	                   const ColumnVector& theta,
	                   mg::Real& o_E) const;


	void computedEdQj(unsigned j,
	                  const mg::Vec3D &m1j,
	                  const mg::Vec3D &m2j,
	                  const ColumnVector& theta,
	                  const mg::Matrix2D& JB,
	                  mg::Real& o_dEQj) const;

	void computeHessian(unsigned j,
	                    const mg::Vec3D &m1j,
	                    const mg::Vec3D &m2j,
	                    const mg::Matrix2D &J,
	                    mg::Real& o_Hjjm1,
	                    mg::Real& o_Hjj,
	                    mg::Real& o_Hjjp1) const;
};
