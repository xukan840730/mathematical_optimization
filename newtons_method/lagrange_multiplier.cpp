#include "../common/common_shared.h"
#include "lagrange_multiplier.h"
#include "newtons_method.h"

struct FuncBodyHelper
{
	ScalarFunc			m_f;
	GradientFunc		m_gf;
	HessianFunc			m_hf;

	ScalarFunc			m_c;
	GradientFunc		m_gc;
	HessianFunc			m_hc;
};

static FuncBodyHelper* g_params0 = nullptr;
static FuncBodyHelper* g_params1 = nullptr;

float PhiFuncBody(const EVector& input)
{
	const FuncBodyHelper* pHelper = static_cast<const FuncBodyHelper*>(g_params0);

	int numParams = input.rows();
	const float lamda = input(numParams - 1);

	EVector gf(numParams - 1);
	EVector gc(numParams - 1);
	pHelper->m_gf(input, &gf);
	pHelper->m_gc(input, &gc);

	EVector ux = gf + gc * lamda;

	float cx = pHelper->m_c(input);

	return ux.squaredNorm() + cx * cx;
}

void gPhiFuncBody(const EVector& input, EVector* output)
{
	xassert(input.rows() == output->rows());

	const FuncBodyHelper* pHelper = static_cast<const FuncBodyHelper*>(g_params1);

	const int numParams = input.rows();
	const float lamda = input(numParams - 1);

	EVector gf(numParams - 1);
	EVector gc(numParams - 1);

	pHelper->m_gf(input, &gf);
	pHelper->m_gc(input, &gc);

	EVector uu = gf + gc * lamda;
	
	EMatrix hf(numParams - 1, numParams - 1);
	EMatrix hc(numParams - 1, numParams - 1);

	pHelper->m_hf(input, &hf);
	pHelper->m_hc(input, &hc);

	EMatrix tt = hf + hc * lamda;

	float cx = pHelper->m_c(input);

	EVector result1 = 2 * tt * uu;
	EVector result2 = 2 * cx * gc;
	EVector result = result1 + result2;

	float dlamba = 2 * (uu.dot(gc));

	// fill results
	xassert(output->rows() == result.rows() + 1);
	for (int ii = 0; ii < numParams - 1; ii++)
		(*output)(ii) = result(ii);

	(*output)(numParams - 1) = dlamba;
}

void LagrangeMultMethod(
	const ScalarFunc F, const GradientFunc gF, const HessianFunc hF,
	const ScalarFunc C, const GradientFunc gC, const HessianFunc hC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, EVector* result)
{
	// valid the constrain function is correct.
	{
		float fC = C(x1);
		xassert(fabsf(fC) < NDI_FLT_EPSILON);
	}

	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	// instead of solving n equations, we optimize the following func:
	// M = (gradient(F(x)) + lamda * (gradient(C(x))))^2 + (C(x))^2

	int numParams = x1.rows();

	ScalarFunc L = PhiFuncBody;
	GradientFunc gL = gPhiFuncBody;

	FuncBodyHelper helper;
	helper.m_f = F;
	helper.m_gf = gF;
	helper.m_hf = hF;
	helper.m_c = C;
	helper.m_gc = gC;
	helper.m_hc = hC;

	g_params0 = &helper;
	g_params1 = &helper;

	NewtonsMethodParams nparams;
	nparams.m_maxIter = params.m_maxIter;
	nparams.m_min = NDI_FLT_EPSILON;

	EVector lx1(numParams + 1);
	for (int ii = 0; ii < numParams; ii++)
		lx1(ii) = x1(ii);
	lx1(numParams) = params.m_lamda1;

	EVector lxstar(numParams + 1);
	QuasiNewtonBFGS(L, gL, nparams, lx1, &lxstar);

	// copy results.
	for (int ii = 0; ii < numParams; ii++)
		(*result)(ii) = lxstar(ii);
}
