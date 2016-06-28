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

void LagrangeMultMethod(
	const ScalarFunc F, const GradientFunc gF, const HessianFunc hF,
	const ScalarFunc C, const GradientFunc gC, const HessianFunc hC,
	const LagrangeMultMethodParams& params,
	const EVector& x1, void* pUserData, void* pReserved, EVector* result)
{
	// valid the constrain function is correct.
	{
		float fC = C(x1, pUserData, pReserved);
		xassert(fabsf(fC) < NDI_FLT_EPSILON);
	}

	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	// instead of solving n equations, we optimize the following func:
	// M = (gradient(F(x)) + lamda * (gradient(C(x))))^2 + (C(x))^2

	int numParams = x1.rows();

	// lagrange function
	auto LFunc = [](const EVector& input, void* pUserData, void* pReserved) -> float {
		xassert(pReserved != nullptr);
		const FuncBodyHelper* pHelper = static_cast<const FuncBodyHelper*>(pReserved);

		int numParams = input.rows();
		const float lamda = input(numParams - 1);

		EVector gf(numParams - 1);
		EVector gc(numParams - 1);
		pHelper->m_gf(input, &gf, pUserData, nullptr);
		pHelper->m_gc(input, &gc, pUserData, nullptr);

		EVector ux = gf + gc * lamda;

		float cx = pHelper->m_c(input, pUserData, nullptr);

		return ux.squaredNorm() + cx * cx;
	};

	// gradiant of lagrange function
	auto gLFunc = [](const EVector& input, EVector* output, void* pUserData, void* pReserved) {
		xassert(pReserved != nullptr);
		xassert(input.rows() == output->rows());

		const FuncBodyHelper* pHelper = static_cast<const FuncBodyHelper*>(pReserved);

		const int numParams = input.rows();
		const float lamda = input(numParams - 1);

		EVector gf(numParams - 1);
		EVector gc(numParams - 1);

		pHelper->m_gf(input, &gf, pUserData, nullptr);
		pHelper->m_gc(input, &gc, pUserData, nullptr);

		EVector uu = gf + gc * lamda;

		EMatrix hf(numParams - 1, numParams - 1);
		EMatrix hc(numParams - 1, numParams - 1);

		pHelper->m_hf(input, &hf);
		pHelper->m_hc(input, &hc);

		EMatrix tt = hf + hc * lamda;

		float cx = pHelper->m_c(input, pUserData, nullptr);

		EVector result1 = 2 * tt * uu;
		EVector result2 = 2 * cx * gc;
		EVector result = result1 + result2;

		float dlamba = 2 * (uu.dot(gc));

		// fill results
		*output = result;
		output->conservativeResize(numParams);
		(*output)(numParams - 1) = dlamba;
	};


	ScalarFunc L = LFunc;
	GradientFunc gL = gLFunc;

	FuncBodyHelper helper;
	helper.m_f = F;
	helper.m_gf = gF;
	helper.m_hf = hF;
	helper.m_c = C;
	helper.m_gc = gC;
	helper.m_hc = hC;

	NewtonsMethodParams nparams;
	nparams.m_maxIter = params.m_maxIter;
	nparams.m_min = NDI_FLT_EPSILON;

	EVector lx1 = x1;
	lx1.conservativeResize(numParams + 1);
	lx1(numParams) = params.m_lamda1;

	EVector lxstar(numParams + 1);
	QuasiNewtonBFGS(L, gL, nparams, lx1, pUserData, &helper, &lxstar);

	// copy results.
	float lamda = lxstar(lxstar.rows() - 1);

	*result = lxstar;
	result->conservativeResize(numParams);
}
