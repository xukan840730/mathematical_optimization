#ifndef _LAGRANGE_MULTIPLIERS_H_
#define _LAGRANGE_MULTIPLIERS_H_

struct LagrangeMultMethodParams
{
public:
	LagrangeMultMethodParams()
		: m_epsilon1(0.001f)
		, m_epsilon2(0.001f)
		, m_lamda1(1.f)
		, m_maxIter(20)
	{}

	float m_epsilon1;	// gradient epsilon to stop the iteration.
	float m_epsilon2;	
	float m_lamda1;		// init guess of lagrange multiplier
	int m_maxIter;		// max number of iterations
};

//void LagrangeMultMethod(
//	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
//	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
//	const LagrangeMultMethodParams& params,
//	const EVector& x1, EVector* result);

// Sequential Quadratic Programming
// handle equality constraint
void SQP1(const CD2Func& objectiveF, const EVector& x0,	const CD2Func& econstrF, const LagrangeMultMethodParams& params, EVector* result);

// handle inequality constraint
void SQP2(const CD2Func& objectiveF, const EVector& x0, const CD2Func& inconstrF, const LagrangeMultMethodParams& params, EVector* result);
// to support any number of inequality constraints
void SQP3(const CD2Func& objectiveF, const EVector& x0, int numInconstr, const CD2Func* inconstrFs, const LagrangeMultMethodParams& params, EVector* result);

#endif
