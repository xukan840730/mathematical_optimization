
#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

struct NewtonsMethodParams
{
	NewtonsMethodParams()
		: m_epsilon(0.001f)
		, m_min(-100.f)
		, m_maxIter(10)
		, m_v(1.f)
	{}

	float m_epsilon;	// gradient epsilon to stop the iteration.
	float m_min;		// safe guard provided by user to stop the iteration.
	int m_maxIter;		// max number of iterations
	float m_v;			// if H matrix is not a positive definite at k iteration, start adding v * I to hessian matrix to become a hessian.
};

void NewtonsMethod(const ScalarFunc F, const GradientFunc g, const HessianFunc H, const NewtonsMethodParams& params,
	const EVector& x1, EVector* result);

// rank one correction.
//void QuasiNewtonSR1(const ScalarF F, const Gradient* g, const NewtonsMethodParams& params,
//	const ScalarVector& x1, ScalarVector* result);

// Davidon Fletcher and Powell method, rank two correction
void QuasiNewtonDFP(const ScalarFunc F, const GradientFunc g, const NewtonsMethodParams& params,
	const EVector& x1, EVector* result);

// BFGS
void QuasiNewtonBFGS(const ScalarFunc F, const GradientFunc g, const NewtonsMethodParams& params,
	const EVector& x1, EVector* result);