#ifndef _LAGRANGE_MULTIPLIERS_H_
#define _LAGRANGE_MULTIPLIERS_H_


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
// SQP3 is not working yet
void SQP3(const CD2Func& objectiveF, const EVector& x0, int numInconstr, const CD2Func* inconstrFs, const LagrangeMultMethodParams& params, EVector* result);

// to support any number of equality constraints
void SQP4(const CD2Func& objectiveF, const EVector& x0, int numEconstr, const CD2Func* econstrFs, const LagrangeMultMethodParams& params, EVector* result);

#endif
