#include "../common/common_shared.h"
#include "lagrange_multiplier.h"
#include "newtons_method.h"

struct FuncBodyHelper
{
	ScalarFunc			m_f;
	GradientFunc		m_gf;
	const EHessian*		m_hf;

	ScalarFunc			m_c;
	GradientFunc		m_gc;
	const EHessian*		m_hc;
} *g_params;

float LagrangeFuncBody(const EVector& input)
{
	const FuncBodyHelper* pHelper = static_cast<const FuncBodyHelper*>(g_params);

	int numParams = input.rows();
	const float lamda = input(numParams - 1);

	EVector gf(numParams - 1);
	EVector gc(numParams - 1);
	pHelper->m_gf(input, &gf);
	pHelper->m_gc(input, &gc);

	EVector ux = gf + gc * lamda;

	float gx = pHelper->m_c(input);

	return ux.squaredNorm() + gx * gx;
}

void LagrangeMultMethod(const ScalarFunc F, const GradientFunc gF, const ScalarFunc C, const GradientFunc gC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, EVector* result)
{
	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	// instead of solving n equations, we optimize the following func:
	// M = (gradient(F(x)) + lamda * (gradient(C(x))))^2 + (C(x))^2

	int numParams = x1.rows();

	ScalarFunc L = LagrangeFuncBody;
	//EGradient gL(numParams + 1);

	FuncBodyHelper helper;
	helper.m_f = F;
	helper.m_gf = gF;
	//helper.m_hf = &HF;
	helper.m_c = C;
	helper.m_gc = gC;
	//helper.m_hc = &HC;

	g_params = &helper;
		
	NewtonsMethodParams nparams;
	nparams.m_maxIter = params.m_maxIter;
	nparams.m_min = NDI_FLT_EPSILON;
	//nparams.extraParams = &helper;


	EVector lx1(numParams + 1);
	for (int ii = 0; ii < numParams; ii++)
		lx1(ii) = x1(ii);
	lx1(numParams) = params.m_lamda1;

	float ftest = L(lx1);
	
	//EVector lxstar(numParams + 1);
	//QuasiNewtonBFGS(L, &gL, nparams, lx1, &lxstar);

	// copy results.
	//for (int ii = 0; ii < numParams; ii++)
	//	(*result)(ii) = lxstar(ii);
}