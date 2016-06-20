
#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

struct NewtonsMethodParams
{
	NewtonsMethodParams()
		: m_epsilon(0.0001f)
		, m_maxIter(10)
		, m_v(1.f)
	{}

	float m_epsilon;	// gradient epsilon to stop the iteration.
	int m_maxIter;		// max number of iterations
	float m_v;			// if H matrix is not a positive definite at k iteration, start adding v * I to hessian matrix to become a hessian.
};

void NewtonsMethod(const ScalarF F, const Gradient* g, const Hessian* H, const NewtonsMethodParams& params,
	const ScalarVector& x1, ScalarVector* result);