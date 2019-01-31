#include "ElasticRod.h"
#include "Utils.h"
#include <dlib/optimization.h>

typedef dlib::matrix<double> Hessian;

//====================================== MinimizationPImpl definition - minimize Energy with respect to twist angle theta =======================================

class ElasticRod::MinimizationPImpl
{
public:
	static void extractThetaVars(const ElasticRod& rod, const ColumnVector& theta, ColumnVector& o_thetaVars);
	static void constructTheta(const ElasticRod& rod, const ColumnVector& thetaVars, ColumnVector& o_theta);

	double minimize(const ElasticRod& rod, ColumnVector& io_theta);

private:
	struct evaluate
	{
		double operator() (const ColumnVector& theta) const;
		const ElasticRod* m_rod = nullptr;
	};

	struct evaluateGradient
	{
		ColumnVector operator() (const ColumnVector& theta) const;
		const ElasticRod* m_rod = nullptr;
	};

	struct evaluateHessian
	{
		Hessian operator() (const ColumnVector& theta) const;
		const ElasticRod* m_rod = nullptr;
	};

	evaluate m_evaluate;
	evaluateGradient m_evaluateGradient;
	evaluateHessian m_evaluateHessian;
};


//====================================== ElasticRod implementation =======================================

ElasticRod::ElasticRod():
    m_minimization( std::make_unique<MinimizationPImpl>() )
{ }

ElasticRod::~ElasticRod() = default;

ElasticRod::ElasticRod( ElasticRod&& other ) = default;

ElasticRod& ElasticRod::operator=( ElasticRod&& other ) = default;

void ElasticRod::initialize( const ElasticRodParams& params,
                             const std::vector<mg::Vec3D>& restpos,
                             const mg::Vec3D& u0,
                             const std::vector<mg::Vec3D>& pos,
                             const std::vector<mg::Vec3D>& vel,
                             const std::vector<mg::Real>& mass,
                             const ColumnVector& theta,
                             const std::set<unsigned> &isClamped )
{
	assert( pos.size() > 2 );
	assert( pos.size() == restpos.size() );
	assert( pos.size() == vel.size() );
	assert( pos.size() == mass.size() );
	assert( static_cast<unsigned>( theta.size() ) == pos.size() - 1 );

	m_params = &params;
	m_u0 = u0;
	m_ppos = pos;
	m_pvel = vel;
	m_pmass = mass;
	m_theta = theta;
	m_isClamped = isClamped;

	m_edges.resize(m_ppos.size() - 1);
	m_kb.resize(m_edges.size());
	m_m1.resize(m_edges.size());
	m_m2.resize(m_edges.size());

	m_restEdgeL.resize(m_edges.size());
	m_restRegionL.resize(m_edges.size());
	m_restWprev.resize(m_edges.size());
	m_restWnext.resize(m_edges.size());

//    compute edges & lengths for rest shape
	computeEdges(restpos, m_edges);
	computeLengths(m_edges, m_restEdgeL, m_restRegionL);
//    compute kb and material frame (Bishop frame) for rest shape
	computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);
//    precompute material curvature for rest shape
	computeMaterialCurvature(m_kb, m_m1, m_m2, m_restWprev, m_restWnext);

//    initialize current state
	updateCurrentState();

#ifdef DBUGG
//    store elasic force for debug purpose
	m_elasticForce.resize(m_ppos.size(), mg::Vec3D(0, 0, 0));
	accumulateInternalElasticForces(m_elasticForce);
#endif

}

void ElasticRod::getState(ElasticRodState& o_state) const
{
	o_state.m_ppos = m_ppos;
	o_state.m_pvel = m_pvel;
	o_state.m_u0 = m_u0;
}

void ElasticRod::setState(const ElasticRodState& state)
{
	m_ppos = state.m_ppos;
	m_pvel = state.m_pvel;
	m_u0 = state.m_u0;

	computeEdges(m_ppos, m_edges);
	updateCurrentState();
}

void ElasticRod::computeEdges(const std::vector<mg::Vec3D>& vertices,
                              std::vector<mg::Vec3D>& o_edges) const
{
	for (unsigned i = 0; i < vertices.size() - 1; ++i)
	{
		o_edges[i] = vertices[i + 1] - vertices[i];
		assert( std::isfinite(o_edges[i].length_squared()) );
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

	for (unsigned i = 1; i < edges.size(); ++i)
	{
///     NOTE: as metioned in the paper the following formula produces ||kb|| = 2 * tan(phi / 2),
///     where phi is the angle of rotation defined by the 2 edges
///     This holds only if ||edges[i]|| = ||rest_edges[i]||,
///     otherwise errors occur in consequent frame rotations and force calculations
///     in particular this assumption can be violated if the constraint enforcement DOES NOT GUARANTEE inextensibility,
///     which is the case with PBD used below
///     note that tan maps [-pi/2, pi/2] to [-inf, inf], so interestingly for arbitrary real number an arbitrary bounded angle can be extracted
///
		o_kb[i] = 2 * mg::cross( edges[i - 1], edges[i] ) /
		        (m_restEdgeL[i - 1] * m_restEdgeL[i] + mg::dot( edges[i - 1], edges[i]) );
	}
}

void ElasticRod::extractSinAndCos(const double &magnitude,
                                  double &o_sinPhi, double &o_cosPhi) const
{
	o_cosPhi = mg::sqrt_safe(4.0 / (4.0 + magnitude));
	o_sinPhi = mg::sqrt_safe(magnitude / (4.0 + magnitude));
}

void ElasticRod::computeBishopFrame(const mg::Vec3D& u0,
                        const std::vector<mg::Vec3D>& edges,
                        std::vector<mg::Vec3D>& o_kb,
                        std::vector<mg::Vec3D>& o_u,
                        std::vector<mg::Vec3D>& o_v) const
{
	computeKB(edges, o_kb);

	o_u[0] = u0;
	o_v[0] = mg::cross(edges[0], o_u[0]);
	o_v[0].normalize();

	double magnitude;
	double sinPhi, cosPhi;
//    compute Bishop frame for current configuration by parallel transporting u0
	for (unsigned i = 1; i < edges.size(); ++i)
	{
//        here sinPhi and cosPhi are derived from the length of kb ||kb|| = 2 * tan( phi/2 )
		magnitude = mg::dot(o_kb[i], o_kb[i]);
		extractSinAndCos(magnitude, sinPhi, cosPhi);
		assert( cosPhi >= 0 && cosPhi <= 1 );

		if ( (1 - cosPhi) < mg::ERR )
		{
			o_u[i] = o_u[i - 1];
			o_v[i] = o_v[i - 1];
			continue;
		}

//        rotate frame u axis around kb
		mg::Quaternion q(cosPhi, sinPhi * mg::normalize(o_kb[i]));
		mg::Quaternion p(0, o_u[i - 1]);
		p = q * p * mg::conjugate(q);

		o_u[i].set(p[1], p[2], p[3]);
		o_u[i].normalize();
		o_v[i] = mg::cross(edges[i], o_u[i]);
		o_v[i].normalize();
	}
}

void ElasticRod::parallelTransportFrame(const mg::Vec3D& e0, const mg::Vec3D& e1,
                                        mg::Vec3D& io_u) const
{
	mg::Vec3D axis = 2 * mg::cross(e0, e1) /
	        (e0.length() * e1.length() + mg::dot(e0, e1));

//    here sinPhi and cosPhi are derived from the length of axis ||axis|| = 2 * tan( phi/2 )
	double sinPhi, cosPhi;
	double magnitude = mg::dot(axis, axis);
	extractSinAndCos(magnitude, sinPhi, cosPhi);
	assert( cosPhi >= 0 && cosPhi <= 1 );

	if ( (1 - cosPhi) < mg::ERR )
	{
		io_u = mg::cross(e1, io_u);
		io_u = mg::normalize( mg::cross(io_u, e1) );
		return;
	}
	mg::Quaternion q(cosPhi, sinPhi * mg::normalize(axis));
	mg::Quaternion p(0, io_u);
	p = q * p * mg::conjugate(q);

	io_u.set(p[1], p[2], p[3]);
	io_u.normalize();
}

void ElasticRod::computeMaterialFrame(const ColumnVector &theta,
                                      std::vector<mg::Vec3D>& io_m1,
                                      std::vector<mg::Vec3D>& io_m2) const
{
	mg::Real sinQ, cosQ;
	mg::Vec3D m1, m2;
	for (unsigned i = 0; i < io_m1.size(); ++i)
	{
		cosQ = std::cos(theta(i));
		sinQ = std::sqrt(1 - cosQ * cosQ);

		m1 = cosQ * io_m1[i] + sinQ * io_m2[i];
		m2 = -sinQ * io_m1[i] + cosQ * io_m2[i];

		io_m1[i] = m1;
		io_m2[i] = m2;
	}
}


void ElasticRod::applyInternalConstraintsIteration()
{
	mg::Vec3D e;
	mg::Real l, l1, l2;
	for (unsigned i = 0; i < m_ppos.size() - 1; ++i)
	{
		bool clamped_i = m_isClamped.count(i) > 0;
		bool clamped_i1 = m_isClamped.count(i + 1) > 0;

		e = m_ppos[i + 1] - m_ppos[i];
//            approximate e.length() with first order accurate Taylor expansion of square root function in the neightbourhood of (restLength^2)
		l = 1 - 2 * m_restEdgeL[i] * m_restEdgeL[i] / (m_restEdgeL[i] * m_restEdgeL[i] +  mg::dot(e, e));

		if (clamped_i)
		{
			l1 = 0;
			l2 = -l;
		}
		else if (clamped_i1)
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

void ElasticRod::accumulateInternalElasticForces(std::vector<mg::Vec3D>& o_forces)
{
	std::vector<mg::Matrix3D> minusGKB(m_edges.size());
	std::vector<mg::Matrix3D> plusGKB(m_edges.size());
	std::vector<mg::Matrix3D> eqGKB(m_edges.size());
	computeGradientKB(m_kb, m_edges, minusGKB, plusGKB, eqGKB);

	std::vector<mg::Vec3D> minusGH(m_kb.size());
	std::vector<mg::Vec3D> plusGH(m_kb.size());
	std::vector<mg::Vec3D> eqGH(m_kb.size());
	computeGradientHolonomyTerms(m_kb, minusGH, plusGH, eqGH);

	mg::Vec2D wkj;
	mg::Matrix2D J;
	mg::matrix_rotation_2D(J, mg::Constants::pi_over_2());

//    compute dE/dQn
	unsigned n = m_edges.size() - 1;
	mg::Real dEdQn;
	computedEdQj(n, m_m1[n], m_m2[n], m_theta, J * m_params->m_B, dEdQn);

	mg::Matrix23D GW;
	mg::Vec3D GH, term;
	for (unsigned i = 0; i < m_ppos.size(); ++i)
	{
		if (m_isClamped.count(i))
		{
			continue;
		}

		for (unsigned k = std::max((int)i - 1, 1); k < m_edges.size(); ++k)
		{
			computeW(m_kb[k], m_m1[k - 1], m_m2[k - 1], wkj);
			computeGradientCurvature(i, k, k - 1,
			                         minusGKB, plusGKB, eqGKB,
			                         minusGH, plusGH, eqGH,
			                         wkj,
			                         J,
			                         GW);

//            o_forces[i] -= (mg::transpose(GW) * m_params->m_B * (wkj - m_restWprev[k])) / m_restRegionL[k];
//            assert( std::isfinite(o_forces[i].length_squared()) );
			term = (mg::transpose(GW) * m_params->m_B * (wkj - m_restWprev[k]));

			computeW(m_kb[k], m_m1[k], m_m2[k], wkj);
			computeGradientCurvature(i, k, k,
			                         minusGKB, plusGKB, eqGKB,
			                         minusGH, plusGH, eqGH,
			                         wkj,
			                         J,
			                         GW);
//            o_forces[i] -= (mg::transpose(GW) * m_params->m_B * (wkj - m_restWnext[k])) / m_restRegionL[k];
//            assert( std::isfinite(o_forces[i].length_squared()) );
			term += (mg::transpose(GW) * m_params->m_B * (wkj - m_restWnext[k]));
			o_forces[i] -= term / m_restRegionL[k];
			assert( std::isfinite(o_forces[i].length_squared()) );
		}

//    need to add  dE/dQn * gradient holonomy(GH) if we have clamped ends since not all twist angles minimize the energy
		if (m_isClamped.size())
		{
			computeGradientHolonomy(i, n, minusGH, plusGH, eqGH, GH);
			o_forces[i] += dEdQn * GH;
			assert( std::isfinite(o_forces[i].length_squared()) );
		}

//        need to limit the force otherwise when e[i] ~ -e[i - 1] ||kb|| goes to infinity => force goes to infinity
		if (o_forces[i].length_squared() > m_params->m_maxElasticForce * m_params->m_maxElasticForce)
		{
			o_forces[i].normalize();
			o_forces[i] *= m_params->m_maxElasticForce;
		}
	}

#ifdef DBUGG
	m_elasticForce = o_forces;
#endif
}

void ElasticRod::computeGradientKB(const std::vector<mg::Vec3D> &kb,
                                   const std::vector<mg::Vec3D> &edges,
                                   std::vector<mg::Matrix3D>& o_minusGKB,
                                   std::vector<mg::Matrix3D>& o_plusGKB,
                                   std::vector<mg::Matrix3D>& o_eqGKB) const
{
// Compute skew-symmetric matrix 3x3 [e], such that [e] * x = cross( e, x )
	std::vector<mg::Matrix3D> edgeMatrix(edges.size());
	for (unsigned i = 0; i < edges.size(); ++i)
	{
		mg::matrix_skew_symmetric(edgeMatrix[i], edges[i]);
	}

	o_minusGKB[0].zero();
	o_plusGKB[0].zero();
	o_eqGKB[0].zero();

	mg::Real scalarFactor;
	for (unsigned i = 1; i < edges.size(); ++i)
	{
		scalarFactor = (m_restEdgeL[i - 1] * m_restEdgeL[i] + mg::dot( edges[i - 1], edges[i] ));
		assert( std::isfinite(scalarFactor) && fabs(scalarFactor) > 0 );

		o_minusGKB[i] = (2.0 * edgeMatrix[i - 1] + mg::outer( kb[i], edges[i - 1] )) / scalarFactor;
		o_plusGKB[i] = (2.0 * edgeMatrix[i] - mg::outer( kb[i], edges[i] )) / scalarFactor;
		o_eqGKB[i] = -(o_plusGKB[i] + o_minusGKB[i]);

//        o_minusGKB[i] = (2.0 * edgeMatrix[i] + mg::outer( kb[i], edges[i] )) / scalarFactor;
//        o_plusGKB[i] = (2.0 * edgeMatrix[i - 1] - mg::outer( kb[i], edges[i - 1] )) / scalarFactor;
//        o_eqGKB[i] = -(o_minusGKB[i] + o_plusGKB[i]);
	}
}

void ElasticRod::computeGradientHolonomyTerms(const std::vector<mg::Vec3D> &kb,
                                       std::vector<mg::Vec3D>& o_minusGH,
                                       std::vector<mg::Vec3D>& o_plusGH,
                                       std::vector<mg::Vec3D>& o_eqGH) const
{
	o_minusGH[0].zero();
	o_plusGH[0].zero();
	o_eqGH[0].zero();

	for (unsigned i = 1; i < kb.size(); ++i)
	{
		o_minusGH[i] = 0.5 * kb[i] / m_restEdgeL[i - 1];
		o_plusGH[i]  = -0.5 * kb[i] / m_restEdgeL[i];
		o_eqGH[i] = -(o_minusGH[i] + o_plusGH[i]);

		assert( std::isfinite( o_minusGH[i].length_squared() ) );
		assert( std::isfinite( o_plusGH[i].length_squared() ) );
		assert( std::isfinite( o_eqGH[i].length_squared() ) );
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
                                        const mg::Matrix2D &J,
                                        mg::Matrix23D &o_GW) const
{
	assert( k >= (i - 1) && (j == k || j == (k - 1)) && j < m_m1.size() );

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

		if (k == (i - 1))
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
	computeGradientHolonomy(i, j, minusGH, plusGH, eqGH, GH);
	o_GW -= J * mg::outer(wkj, GH);
}

void ElasticRod::computeGradientHolonomy(unsigned i , unsigned j,
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
	if (j >= i && i < eqGH.size())
	{
		o_GH += eqGH[i];
	}
	if (j >= (i + 1) && (i + 1) < minusGH.size())
	{
		o_GH += minusGH[i + 1];
	}
}

void ElasticRod::computeEnergy(const std::vector<mg::Vec3D>& m1,
                               const std::vector<mg::Vec3D>& m2,
                               const ColumnVector &theta,
                               mg::Real &o_E) const
{
	o_E = 0.0;
	mg::Vec2D wij;
	mg::Real mi;
	for (unsigned i = 1; i < m_edges.size(); ++i)
	{
//        bend energy term
		computeW(m_kb[i], m1[i - 1], m2[i - 1], wij);
		o_E += mg::dot(wij - m_restWprev[i], m_params->m_B * (wij - m_restWprev[i])) * 0.5 / m_restRegionL[i];

		computeW(m_kb[i], m1[i], m2[i], wij);
		o_E += mg::dot(wij - m_restWnext[i], m_params->m_B * (wij - m_restWnext[i])) * 0.5 / m_restRegionL[i];

//        twist energy term
		mi = (theta(i) - theta(i - 1));
		o_E += m_params->m_beta * mi * mi / m_restRegionL[i];
	}
}

void ElasticRod::computedEdQj(unsigned j,
                              const mg::Vec3D& m1j,
                              const mg::Vec3D& m2j,
                              const ColumnVector &theta,
                              const mg::Matrix2D &JB,
                              mg::Real &o_dEQj) const
{
	o_dEQj = 0.0;
	mg::Vec2D wij;
	mg::Real term;
//    compute first term dWj/dQj + 2 * beta * mj / lj
	if (j > 0)
	{
		computeW(m_kb[j], m1j, m2j, wij);

		term = mg::dot(wij, JB * (wij - m_restWnext[j]));
		term += 2 * m_params->m_beta * (theta(j) - theta(j - 1));
		term /= m_restRegionL[j];

		o_dEQj += term;
	}
//    compute second term dWj+1/dQj - 2 * beta * mj+1 / lj+1
	if (j < m_edges.size() - 1)
	{
		computeW(m_kb[j + 1], m1j, m2j, wij);

		term = mg::dot(wij, JB * (wij - m_restWprev[j + 1]));
		term -= 2 * m_params->m_beta * (theta(j + 1) - theta(j));
		term /= m_restRegionL[j + 1];

		o_dEQj += term;
	}
}

void ElasticRod::computeHessian(unsigned j,
                                const mg::Vec3D& m1j,
                                const mg::Vec3D& m2j,
                                const mg::Matrix2D& J,
                                mg::Real &o_Hjjm1,
                                mg::Real &o_Hjj,
                                mg::Real &o_Hjjp1) const
{
	o_Hjjm1 = o_Hjj = o_Hjjp1 = 0;

	mg::Vec2D wij;
	double hjj;
	if (j > 0)
	{
		o_Hjjm1 = -2 * m_params->m_beta / m_restRegionL[j];

		computeW(m_kb[j], m1j, m2j, wij);

		hjj = 2 * m_params->m_beta;
		hjj += mg::dot( wij, mg::transpose(J) * m_params->m_B * J * wij );
		hjj -= mg::dot( wij, m_params->m_B * (wij - m_restWnext[j]) );
		hjj /= m_restRegionL[j];

		o_Hjj = hjj;
	}
	if (j < m_edges.size() - 1)
	{
		o_Hjjp1 = -2 * m_params->m_beta / m_restRegionL[j + 1];

		computeW(m_kb[j + 1], m1j, m2j, wij);
		hjj = 2 * m_params->m_beta;
		hjj += mg::dot( wij, mg::transpose(J) * m_params->m_B * J * wij );
		hjj -= mg::dot( wij, m_params->m_B * (wij - m_restWprev[j + 1]) );
		hjj /= m_restRegionL[j + 1];

		o_Hjj += hjj;
	}
}

void ElasticRod::updateCurrentState()
{
//    parallel transport first frame in time
	mg::Vec3D e0 = m_edges[0];
	computeEdges(m_ppos, m_edges);
	mg::Vec3D e1 = m_edges[0];

	parallelTransportFrame(e0, e1, m_u0);
	computeBishopFrame(m_u0, m_edges, m_kb, m_m1, m_m2);

	double minE = 0;
	if (m_params->m_strategy != ElasticRodParams::MinimizationStrategy::NONE)
	{
		minE = m_minimization->minimize(*this, m_theta);
	}
	computeMaterialFrame(m_theta, m_m1, m_m2);

#ifdef DBUGG
	mg::Real E;
	computeEnergy(m_m1, m_m2, m_theta, E);
	std::cout<< "Theta:" << m_theta << "\n";
	std::cout<< "Total Energy: " << E << " MIN Energy: " << minE << std::endl;
#endif
}




//====================================== MinimizationPImpl implementation =======================================

void ElasticRod::MinimizationPImpl::extractThetaVars(const ElasticRod& rod, const ColumnVector& theta, ColumnVector& o_thetaVars )
{
	o_thetaVars.set_size( theta.size() - rod.m_isClamped.size() );
	auto j = 0l;
	for ( auto i = 0l; i < theta.size(); ++i )
	{
		if ( rod.m_isClamped.count( i ) || rod.m_isClamped.count( i + 1 ) )
		{
			continue;
		}
		o_thetaVars( j ) = theta( i );
		++j;
	}
}

void ElasticRod::MinimizationPImpl::constructTheta( const ElasticRod& rod, const ColumnVector& thetaVars, ColumnVector& o_theta )
{
	auto j = 0l;
	for ( auto i = 0l; i < o_theta.size(); ++i )
	{
		if (rod.m_isClamped.count( i ) || rod.m_isClamped.count( i + 1 ) )
		{
			continue;
		}
		o_theta( i ) = thetaVars( j );
		++j;
	}
}

double ElasticRod::MinimizationPImpl::minimize(const ElasticRod& rod, ColumnVector& io_theta)
{
	m_evaluate.m_rod = &rod;
	m_evaluateGradient.m_rod = &rod;
	m_evaluateHessian.m_rod = &rod;

	ColumnVector thetaVars;
	extractThetaVars(rod, io_theta, thetaVars);

	double minE = -1;
	switch (rod.m_params->m_strategy) {
	case ElasticRodParams::MinimizationStrategy::BFGS_NUMERIC:
		minE = dlib::find_min(dlib::bfgs_search_strategy(),
		                     dlib::objective_delta_stop_strategy(rod.m_params->m_tolerance, rod.m_params->m_maxIter),
		                     m_evaluate,
		                     dlib::derivative(m_evaluate),
		                     thetaVars,
		                     0.0);
		break;
	case ElasticRodParams::MinimizationStrategy::BFGS:
		minE = dlib::find_min(dlib::bfgs_search_strategy(),
		                     dlib::objective_delta_stop_strategy(rod.m_params->m_tolerance, rod.m_params->m_maxIter),
		                     m_evaluate,
		                     m_evaluateGradient,
		                     thetaVars,
		                     0.0);
		break;
	case ElasticRodParams::MinimizationStrategy::NEWTON:
		minE = dlib::find_min(dlib::newton_search_strategy(m_evaluateHessian),
		                     dlib::objective_delta_stop_strategy(rod.m_params->m_tolerance, rod.m_params->m_maxIter),
		                     m_evaluate,
		                     m_evaluateGradient,
		                     thetaVars,
		                     0.0);
		break;
	default:
		break;
	}

	constructTheta(rod, thetaVars, io_theta);

	return minE;
}

double ElasticRod::MinimizationPImpl::evaluate::operator ()(const ColumnVector& theta) const
{
	assert( m_rod != nullptr );

//    TODO FIX: can I avoid copying ?????
//    keeping m1, m2 as members and copying the data gives ~ 1ms performance benefit
	std::vector<mg::Vec3D> m1 = m_rod->m_m1;
	std::vector<mg::Vec3D> m2 = m_rod->m_m2;
	ColumnVector theta_full = m_rod->m_theta;

	constructTheta( *m_rod, theta, theta_full );
	m_rod->computeMaterialFrame( theta_full, m1, m2 );

	mg::Real E;
	m_rod->computeEnergy( m1, m2, theta_full, E );
	return static_cast<double>( E );
}

ColumnVector ElasticRod::MinimizationPImpl::evaluateGradient::operator ()( const ColumnVector &theta ) const
{
	assert( m_rod != nullptr );

	ColumnVector gradient( theta.size() );
//    TODO FIX: can I avoid copying ?????
	std::vector<mg::Vec3D> m1 = m_rod->m_m1;
	std::vector<mg::Vec3D> m2 = m_rod->m_m2;
	ColumnVector theta_full = m_rod->m_theta;

	constructTheta( *m_rod, theta, theta_full );
	m_rod->computeMaterialFrame( theta_full, m1, m2 );

	mg::Matrix2D JB;
	mg::matrix_rotation_2D( JB, mg::Constants::pi_over_2() );
	JB *= m_rod->m_params->m_B;

	auto j = 0ll;
	for ( auto i = 0ll; i < theta_full.size(); ++i )
	{
		if ( m_rod->m_isClamped.count( i ) || m_rod->m_isClamped.count( i + 1 ) )
		{
			continue;
		}
		mg::Real dEdQj;
		m_rod->computedEdQj( i, m1[i], m2[i], theta_full, JB, dEdQj );
		gradient( j ) = dEdQj;
		++j;
	}

	return gradient;
}

Hessian ElasticRod::MinimizationPImpl::evaluateHessian::operator ()( const ColumnVector& theta ) const
{
	assert( m_rod != nullptr );

	Hessian hessian( theta.size(), theta.size() );
	hessian = dlib::zeros_matrix( hessian );
//    TODO FIX: can I avoid copying ?????
	std::vector<mg::Vec3D> m1 = m_rod->m_m1;
	std::vector<mg::Vec3D> m2 = m_rod->m_m2;
	ColumnVector theta_full = m_rod->m_theta;

	constructTheta( *m_rod, theta, theta_full );
	m_rod->computeMaterialFrame( theta_full, m1, m2 );

	mg::Matrix2D J;
	mg::matrix_rotation_2D( J, mg::Constants::pi_over_2() );
	auto j = 0ll;
	for ( auto i = 0ll; i < m_rod->m_theta.size(); ++i )
	{
		if ( m_rod->m_isClamped.count( i ) || m_rod->m_isClamped.count( i + 1 ) )
		{
			continue;
		}

		mg::Real Hjjm1, Hjj, Hjjp1;
		m_rod->computeHessian( i, m1[i], m2[i], J, Hjjm1, Hjj, Hjjp1 );

		if ( j > 0 )
		{
			hessian( j, j - 1 ) = Hjjm1;
		}
		hessian( j, j ) = Hjj;
		if ( ( j + 1 ) < theta.size() )
		{
			hessian( j, j + 1 ) = Hjjp1;
		}
		++j;
	}

	return hessian;
}









//============================  =======================



//    mg::Real dEdQj;
//        for (unsigned j = std::max((int)(i) - 1, 0); j < m_edges.size(); ++j)
//        {
//            computedEdQj(j, J * m_B, dEdQj);
//            computeGradientHolonomySum(i, n, minusGH, plusGH, eqGH, GH);
//            o_forces[i] += dEdQj * GH;
//        }


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

//        if (m_isClamped.count(i))
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

