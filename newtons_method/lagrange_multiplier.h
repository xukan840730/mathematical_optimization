#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

struct LagrangeMultMethodParams
{
public:
	LagrangeMultMethodParams()
		: m_lamda1(1.f)
		, m_maxIter(20)
	{}

	float m_lamda1;		// init guess of lagrange multiplier
	int m_maxIter;		// max number of iterations
};

// Lagrange Multipliers Method, solve Lagrange function numerically
void LagrangeMultMethod(
	const ScalarFunc F, const GradientFunc gF, const HessianFunc hF,
	const ScalarFunc C, const GradientFunc gC, const HessianFunc hC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, EVector* result);