#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

struct LagrangeMultMethodParams
{
public:
	LagrangeMultMethodParams()
		: m_epsilon1(0.001f)
		, m_epsilon2(0.001f)
		, m_lamda1(0.f)
		, m_maxIter(20)
	{}

	float m_epsilon1;	// gradient epsilon to stop the iteration.
	float m_epsilon2;
	float m_lamda1;		// init guess of lagrange multiplier
	int m_maxIter;		// max number of iterations
};

struct LagrangeMultMethodResult
{
	LagrangeMultMethodResult()
		: m_iter(0)
	{}

	LagrangeMultMethodResult(int iter)
		: m_iter(iter)
	{}

	int m_iter;
};

// Lagrange Multipliers Method, solve Lagrange function numerically
LagrangeMultMethodResult LagrangeMultMethod(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, EVector* result);


// Sequential Quadratic Programming
void SQP1(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, EVector* result);