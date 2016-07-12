#include "../common/common_shared.h"
#include "lagrange_multiplier.h"
#include "../unconstrained_optimization/newtons_method.h"

LagrangeMultMethodResult LagrangeMultMethod(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
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

	const float scale = 1.f;

	// lagrange function
	auto LFunc = [&F, &gF, &C, &gC, scale](const EVector& input) -> float {
		int numParams = input.rows();
		const float lamda = input(numParams - 1);

		EVector gf(numParams - 1);
		EVector gc(numParams - 1);
		gF(input, &gf);
		gC(input, &gc);

		EVector ux = gf + gc * lamda * scale;

		float cx = C(input);

		return ux.squaredNorm() + cx * cx * scale * scale;
	};

	// gradiant of lagrange function
	auto gLFunc = [&F, &gF, &hF, &C, &gC, &hC, scale](const EVector& input, EVector* output) {
		xassert(input.rows() == output->rows());

		const int numParams = input.rows();
		const float lamda = input(numParams - 1);

		EVector gf(numParams - 1);
		EVector gc(numParams - 1);

		gF(input, &gf);
		gC(input, &gc);

		EVector uu = gf + gc * lamda * scale;

		EMatrix hf(numParams - 1, numParams - 1);
		EMatrix hc(numParams - 1, numParams - 1);

		hF(input, &hf);
		hC(input, &hc);

		EMatrix tt = hf + hc * lamda * scale;

		float cx = C(input);

		EVector result1 = 2 * tt * uu;
		EVector result2 = 2 * scale * scale * cx * gc;
		EVector result = result1 + result2;

		float dlamba = 2 * scale * (uu.dot(gc));

		// fill results
		*output = result;
		output->conservativeResize(numParams);
		(*output)(numParams - 1) = dlamba;
	};


	ScalarFunc L = LFunc;
	GradientFunc gL = gLFunc;

	NewtonsMethodParams nparams;
	nparams.m_maxIter = params.m_maxIter;
	nparams.m_min = NDI_FLT_EPSILON;

	EVector lx1 = x1;
	lx1.conservativeResize(numParams + 1);
	lx1(numParams) = params.m_lamda1;

	EVector lxstar(numParams + 1);
	NewtonsMethodResult nres = QuasiNewtonBFGS(L, gL, nparams, lx1, &lxstar);

	// copy results.
	float lamda = lxstar(lxstar.rows() - 1);

	*result = lxstar;
	result->conservativeResize(numParams);
	
	LagrangeMultMethodResult res;
	res.m_iter = nres.m_iter;
	return res;
}

//-----------------------------------------------------------------------------------------------------------//
void SQP1(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
	const LagrangeMultMethodParams& params, const EVector& x1, EVector* result)
{
	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	const int numParams = x1.rows();

	// L func
	auto LFunc = [&F, &C](const EVector& input, const float lambda) -> float {
		return F(input) + C(input) * lambda;
	};

	// gradient of L func respect to x vector
	auto gLFunc = [&gF, &gC, numParams](const EVector& input, const float lambda, EVector* output) {
		EVector res1(numParams);
		EVector res2(numParams);
		gF(input, &res1);
		gC(input, &res2);
		*output = res1 + res2 * lambda;
	};

	// hessian of L func respect to x vector
	auto hLFunc = [&hF, &hC, numParams](const EVector& input, const float lambda, EMatrix* output) {
		EMatrix res1(numParams, numParams);
		EMatrix res2(numParams, numParams);
		hF(input, &res1);
		hC(input, &res2);
		*output = res1 + res2 * lambda;
	};

	float lambdaK = params.m_lamda1;
	EVector xk = x1;

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		// evaluate g(k)
		EVector gk(numParams);
		gF(xk, &gk);

		{
			float norm2 = gk.squaredNorm();
			if (norm2 < params.m_epsilon1 * params.m_epsilon1)
			{
				break;
			}
		}

		EMatrix A;
		{
			EMatrix hLk(numParams, numParams);
			hLFunc(xk, lambdaK, &hLk);
			A = hLk;

			EVector gCk(numParams);
			gC(xk, &gCk);

			A.conservativeResize(numParams + 1, numParams + 1);
			// TODO: replace with Eigen interface
			for (int i = 0; i < numParams; i++)
			{
				A(numParams, i) = A(i, numParams) = gCk(i);
			}
			A(numParams, numParams) = 0.f;
		}

		EVector B;
		{
			gLFunc(xk, lambdaK, &B);
			B.conservativeResize(numParams + 1);

			float ck = C(xk);
			B(numParams) = ck;

			B *= -1.f;
		}

		// solve A * delta = B;
		EMatrix deltaM = A.fullPivLu().solve(B);
		EVector deltaV = deltaM.col(0);

		//xk += deltaV;
		{
			for (int ii = 0; ii < xk.rows(); ii++)
			{
				xk(ii) += deltaV(ii);
			}
		}
		lambdaK += deltaV(numParams);

		if (deltaV.squaredNorm() < params.m_epsilon2 * params.m_epsilon2)
		{
			break;
		}
	}

	*result = xk;
}
