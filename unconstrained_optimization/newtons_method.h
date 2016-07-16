
#ifndef _NEWTONS_METHOD_H_
#define _NEWTONS_METHOD_H_

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

// classical newton's method.
void NewtonsMethod(const CD2Func& objectiveF, const EVector& x1, const NewtonsMethodParams& params, EVector* result);

// BFGS
void QuasiNewtonBFGS(const CD1Func& objectiveF, const EVector& x1, const NewtonsMethodParams& params, EVector* result);

#endif