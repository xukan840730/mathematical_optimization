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
// Lagrange functions:
typedef std::function<float(const EVector& input, float lambda)> LScalarFunc;
typedef std::function<void(const EVector& input, float lambda, EVector* output)> LGradientFunc;
typedef std::function<void(const EVector& input, float lambda, EMatrix* output)> LHessianFunc;

//-----------------------------------------------------------------------------------------------------------//
static void SQPFramework(
	const LScalarFunc& LFunc, const LGradientFunc& gLFunc, const LHessianFunc& hLFunc,
	const ScalarFunc& HFunc, const GradientFunc& gHFunc, //const HessianFunc& hHFunc,
	int maxIter, float lambda0, float epsilon1, float epsilon2, EVector x0, EVector* result)
{
	const int numVars = x0.rows();
	EVector xk = x0;
	float lambdaK = lambda0;

	for (int iter = 1; iter <= maxIter; iter++)
	{
		// evaluate g(k)
		EVector gk(numVars);
		gLFunc(xk, lambdaK, &gk);

		{
			float norm2 = gk.squaredNorm();
			if (norm2 < epsilon1 * epsilon1)
			{
				break;
			}
		}

		EMatrix A;
		{
			EMatrix hLk(numVars, numVars);
			hLFunc(xk, lambdaK, &hLk);
			A = hLk;

			EVector gCk(numVars);
			gHFunc(xk, &gCk);

			A.conservativeResize(numVars + 1, numVars + 1);
			// TODO: replace with Eigen interface
			for (int i = 0; i < numVars; i++)
				A(numVars, i) = A(i, numVars) = gCk(i);

			A(numVars, numVars) = 0.f;
		}

		EVector B;
		{
			gLFunc(xk, lambdaK, &B);
			B.conservativeResize(numVars + 1);

			float ck = HFunc(xk);
			B(numVars) = ck;

			B *= -1.f;
		}

		// solve A * delta = B;
		EMatrix deltaM = A.fullPivLu().solve(B);
		EVector deltaV = deltaM.col(0);

		//xk += deltaV;
		{
			for (int ii = 0; ii < xk.rows(); ii++)
				xk(ii) += deltaV(ii);
		}
		lambdaK += deltaV(numVars);

		if (deltaV.squaredNorm() < epsilon2 * epsilon2)
		{
			break;
		}
	}

	*result = xk;

}

//-----------------------------------------------------------------------------------------------------------//
void SQP1(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
	const LagrangeMultMethodParams& params, const EVector& x0, EVector* result)
{
	const int numOrigParams = x0.rows();

	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	// make initial lamda0 guess.
	// lambdaStar = -[gC(xstar)^t . A . gC(xstar)]^-1 . gC(xstar)^t . A . gF(xstar), where A is nonsingular matrix that is positive definite on null space of gC(xstar)^t
	// let A be I, lambda0 = -[gC(x0)^t * gC(x0)]^-1 . gC(x0)^t . gF(x0)

	float lambda0;
	{
		EVector gC0(numOrigParams);
		gC(x0, &gC0);
		EVector gF0(numOrigParams);
		gF(x0, &gF0);

		float uu = gC0.transpose().dot(gC0);
		lambda0 = -(gC0.transpose().dot(gF0)) / uu;
	}

	// L func
	auto LFunc = [&F, &C](const EVector& input, const float lambda) -> float {
		return F(input) + C(input) * lambda;
	};

	// gradient of L func respect to x vector
	auto gLFunc = [&gF, &gC, numOrigParams](const EVector& input, const float lambda, EVector* output) {
		EVector res1(numOrigParams);
		EVector res2(numOrigParams);
		gF(input, &res1);
		gC(input, &res2);
		*output = res1 + res2 * lambda;
	};

	// hessian of L func respect to x vector
	auto hLFunc = [&hF, &hC, numOrigParams](const EVector& input, const float lambda, EMatrix* output) {
		EMatrix res1(numOrigParams, numOrigParams);
		EMatrix res2(numOrigParams, numOrigParams);
		hF(input, &res1);
		hC(input, &res2);
		*output = res1 + res2 * lambda;
	};

	SQPFramework(LFunc, gLFunc, hLFunc, C, gC, //hC, 
		params.m_maxIter, lambda0, params.m_epsilon1, params.m_epsilon2, x0, result);
}

//-----------------------------------------------------------------------------------------------------------//
void SQP2(
	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
	const LagrangeMultMethodParams& params, const EVector& x0, EVector* result)
{
	// change inequality constraint to equality constraint by introducing slack variable.
	// oinput is the original variables vector
	// ninput is the variables vector with slack variables.

	auto nF = [&F](const EVector& ninput) -> float {
		EVector oinput = ninput;
		oinput.conservativeResize(ninput.rows() - 1);
		return F(oinput);
	};

	auto gnF = [gF](const EVector& ninput, EVector* noutput) {
		int numORows = ninput.rows() - 1;

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EVector ooutput(numORows);
		gF(oinput, &ooutput);
		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1);
		(*noutput)(numORows) = 0.f;
	};

	auto hnF = [&hF](const EVector& ninput, EMatrix* noutput) {
		int numORows = ninput.rows() - 1;
		
		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EMatrix ooutput(numORows, numORows);
		hF(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1, numORows + 1);

		for (int ii = 0; ii < numORows; ii++)
			(*noutput)(ii, numORows) = (*noutput)(numORows, ii) = 0.f;

		(*noutput)(numORows, numORows) = 0.f;
	};

	// h(x) = C(x) + slk * slk == 0
	auto HFunc = [&C](const EVector& ninput) -> float {
		int numORows = ninput.rows() - 1;
		float slackVar = ninput(numORows);

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);
		return C(oinput) + slackVar * slackVar;
	};

	// gh(x) = (gC(x), 2 * slk)
	auto gHFunc = [&gC](const EVector& ninput, EVector* noutput) {
		int numORows = ninput.rows() - 1;
		float slackVar = ninput(numORows);

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);
		
		EVector ooutput(numORows);
		gC(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1);
		(*noutput)(numORows) = 2.f * slackVar;
	};

	// hh(x) = ...
	auto hHFunc = [&hC](const EVector& ninput, EMatrix* noutput) {
		int numORows = ninput.rows() - 1;

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EMatrix ooutput(numORows, numORows);
		hC(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1, numORows + 1);

		for (int ii = 0; ii < numORows; ii++)
			(*noutput)(ii, numORows) = (*noutput)(numORows, ii) = 0.f;

		(*noutput)(numORows, numORows) = 2.f;
	};

	// new L func with converted equality constraint
	auto LFunc = [&nF, &HFunc](const EVector& ninput, const float lambda) -> float {
		int numORows = ninput.rows() - 1;

		EVector oinput = ninput;
		oinput.conservativeResize(numORows - 1);

		return nF(oinput) + HFunc(ninput) * lambda;
	};

	// gradient of L func respect to x vector
	auto gLFunc = [&gnF, &gHFunc](const EVector& ninput, const float lambda, EVector* output) {
		int numNVars = ninput.rows();
		EVector res1(numNVars);
		EVector res2(numNVars);
		gnF(ninput, &res1);
		gHFunc(ninput, &res2);
		*output = res1 + res2 * lambda;
	};

	// hessian of L func respect to x vector
	auto hLFunc = [&hnF, &hHFunc](const EVector& ninput, const float lambda, EMatrix* output) {
		int numNVars = ninput.rows();
		EMatrix res1(numNVars, numNVars);
		EMatrix res2(numNVars, numNVars);
		hnF(ninput, &res1);
		hHFunc(ninput, &res2);
		*output = res1 + res2 * lambda;
	};

	const int numOVars = x0.rows();
	const int numNVars = numOVars + 1;

	float lambda0 = params.m_lamda1;
	EVector nx0 = x0;
	nx0.conservativeResize(numOVars + 1);
	nx0(numOVars) = 0.f;
	//{
	//	float initC = _C(x0);
	//	ASSERTF(initC <= 0.f, ("initial guess x0 doesn't satisfy constraint!"));
	//	if (initC >= 0.f)
	//	{
	//		return;
	//	}
	//	xk(numOVars) = sqrtf(-initC);
	//}

	EVector nresult(numNVars);
	SQPFramework(
		LFunc, gLFunc, hLFunc, 
		HFunc, 
		gHFunc, 
		//hHFunc,
		params.m_maxIter, 
		lambda0, 
		params.m_epsilon1, 
		params.m_epsilon2, 
		nx0, 
		&nresult);

	*result = nresult;
	result->conservativeResize(numOVars);
}