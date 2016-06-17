
#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

struct NewtonsMethodParams
{
	NewtonsMethodParams()
		: m_epsilon(0.0001f)
		, m_maxIter(10)
	{}

	float m_epsilon;	// gradient epsilon to stop the iteration.
	int m_maxIter;		// max number of iterations
};

void NewtonsMethod(const ScalarF F, const Gradient* g, const Hessian* H, const NewtonsMethodParams& params,
	const ScalarVector& initGuess, ScalarVector* result);