#include "../common/common_shared.h"
#include "lagrange_multiplier.h"
//#include "ndlib/math/newtons_method.h"

//void LagrangeMultMethod(
//	const ScalarFunc& F, const GradientFunc& gF, const HessianFunc& hF,
//	const ScalarFunc& C, const GradientFunc& gC, const HessianFunc& hC,
//	const LagrangeMultMethodParams& params,
//	const EVector& x1, EVector* result)
//{
//	// valid the constrain function is correct.
//	{
//		float fC = C(x1);
//		ALWAYS_ASSERT(fabsf(fC) < NDI_FLT_EPSILON);
//	}
//
//	// L = F(x) + lamda * C(x)
//	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0
//
//	// instead of solving n equations, we optimize the following func:
//	// M = (gradient(F(x)) + lamda * (gradient(C(x))))^2 + (C(x))^2
//
//	int numParams = x1.rows();
//
//	// lagrange function
//	auto LFunc = [&F, &gF, &C, &gC](const EVector& input) -> float {
//		int numParams0 = input.rows();
//		const float lamda = input(numParams0 - 1);
//
//		EVector gf(numParams0 - 1);
//		EVector gc(numParams0 - 1);
//		gF(input, &gf);
//		gC(input, &gc);
//
//		EVector ux = gf + gc * lamda;
//
//		float cx = C(input);
//
//		return ux.squaredNorm() + cx * cx;
//	};
//
//	// gradiant of lagrange function
//	auto gLFunc = [&F, &gF, &hF, &C, &gC, &hC](const EVector& input, EVector* output) {
//		ALWAYS_ASSERT(input.rows() == output->rows());
//
//		const int numParams0 = input.rows();
//		const float lamda = input(numParams0 - 1);
//
//		EVector gf(numParams0 - 1);
//		EVector gc(numParams0 - 1);
//
//		gF(input, &gf);
//		gC(input, &gc);
//
//		EVector uu = gf + gc * lamda;
//
//		EMatrix hf(numParams0 - 1, numParams0 - 1);
//		EMatrix hc(numParams0 - 1, numParams0 - 1);
//
//		hF(input, &hf);
//		hC(input, &hc);
//
//		EMatrix tt = hf + hc * lamda;
//
//		float cx = C(input);
//
//		EVector result1 = 2 * tt * uu;
//		EVector result2 = 2 * cx * gc;
//		EVector finalRes = result1 + result2;
//
//		float dlamba = 2 * (uu.dot(gc));
//
//		// fill results
//		*output = finalRes;
//		output->conservativeResize(numParams0);
//		(*output)(numParams0 - 1) = dlamba;
//	};
//
//	CD1Func L(LFunc, gLFunc);
//
//	NewtonsMethodParams nparams;
//	nparams.m_maxIter = params.m_maxIter;
//	nparams.m_min = NDI_FLT_EPSILON;
//
//	EVector lx1 = x1;
//	lx1.conservativeResize(numParams + 1);
//	lx1(numParams) = params.m_lamda1;
//
//	EVector lxstar(numParams + 1);
//	QuasiNewtonBFGS(L, lx1, nparams, &lxstar);
//
//	// copy results.
//	float lamda = lxstar(lxstar.rows() - 1);
//
//	*result = lxstar;
//	result->conservativeResize(numParams);
//}

//-----------------------------------------------------------------------------------------------------------//
// helper function for eigen vectors.
//-----------------------------------------------------------------------------------------------------------//
static void ChangeEVector(EVector* inout, int numRows)
{
	int numORows = inout->rows();
	if (numORows < numRows)
	{
		inout->conservativeResize(numRows);
		for (int i = numORows; i < numRows; i++)
			(*inout)(i) = 0.f;
	}
	else if (numORows > numRows)
	{
		inout->conservativeResize(numRows);
	}
}

static void ChangeEMatrix(EMatrix* inout, int numRows)
{
	int numORows = inout->rows();
	if (numORows < numRows)
	{
		inout->conservativeResize(numRows, numRows);

		for (int ii = numORows; ii < numRows; ii++)
			for (int jj = numRows; jj < numRows; jj++)
				(*inout)(ii, jj) = (*inout)(jj, ii) = 0.f;

		(*inout)(numORows, numORows) = 0.f;
	}
	else
	{
		inout->conservativeResize(numRows, numRows);
	}
}

static EVector EVectorTimesLambda(const EVector& input, const EVector& lambdas)
{
	ASSERT(input.rows() == lambdas.rows());
	EVector result(input.rows());

	for (int i = 0; i < input.rows(); i++)
		result(i) = input(i) * lambdas(i);

	return result;
}

static EMatrix EMatrixTimesLambda(const EMatrix& input, const EVector& lambdas)
{
	ASSERT(input.rows() == lambdas.rows());
	EMatrix result(input.rows(), input.rows());

	for (int i = 0; i < input.rows(); i++)
		result.row(i) *= lambdas(i);

	return result;
}

//-----------------------------------------------------------------------------------------------------------//
// Lagrange functions:
typedef std::function<float(const EVector& input, const EVector& lambdas)> LScalarFunc;
typedef std::function<void(const EVector& input, const EVector& lambdas, EVector* output)> LGradientFunc;
typedef std::function<void(const EVector& input, const EVector& lambdas, EMatrix* output)> LHessianFunc;

// 1 time continous differentiable lagrange function
struct L1Func {
	const LScalarFunc& f;
	const LGradientFunc& g;

	L1Func(const LScalarFunc& _f, const LGradientFunc& _g) : f(_f), g(_g) {}
};

// 2 times continous differentiable lagrange function
struct L2Func {
	const LScalarFunc& f;
	const LGradientFunc& g;
	const LHessianFunc& h;

	L2Func(const LScalarFunc& _f, const LGradientFunc& _g, const LHessianFunc& _h) : f(_f) , g(_g) , h(_h) {}
};

//-----------------------------------------------------------------------------------------------------------//
static void SQPFramework(const L2Func& LF, const EVector& x0, const EVector& lambda0,
	int numEConstr, const CD1Func* econstrFs, int maxIter, float epsilon1, float epsilon2, EVector* result)
{
	const int numOVars = x0.rows();
	EVector xk = x0;
	EVector lambdaK = lambda0;
	ASSERT(numEConstr == lambdaK.rows());

	LF.f(xk, lambdaK);

	for (int iter = 1; iter <= maxIter; iter++)
	{
		// evaluate g(k)
		EVector gk(numOVars);
		LF.g(xk, lambdaK, &gk);

		{
			float norm2 = gk.squaredNorm();
			if (norm2 < epsilon1 * epsilon1)
			{
				break;
			}
		}

		EMatrix A;
		{
			EMatrix hLk(numOVars, numOVars);
			LF.h(xk, lambdaK, &hLk);
			A = hLk;

			A.conservativeResize(numOVars + numEConstr, numOVars + numEConstr);
			// TODO: replace with Eigen interface
			for (int iConstr = 0; iConstr < numEConstr; iConstr++)
			{
				EVector gCk(numOVars);
				econstrFs[iConstr].g(xk, &gCk);

				for (int jj = 0; jj < numOVars; jj++)
					A(numOVars + iConstr, jj) = A(jj, numOVars + iConstr) = gCk(jj);
			}

			for (int iConstr = 0; iConstr < numEConstr; iConstr++)
				A(numOVars + iConstr, numOVars + iConstr) = 0.f;
		}

		EVector B;
		{
			LF.g(xk, lambdaK, &B);
			B.conservativeResize(numOVars + numEConstr);

			for (int ii = 0; ii < numEConstr; ii++)
			{
				float ck = econstrFs[ii].f(xk);
				B(numOVars + ii) = ck;
			}

			B *= -1.f;
		}

		// solve A * delta = B;
		EMatrix deltaM = A.fullPivLu().solve(B);
		EVector deltaV = deltaM.col(0);

		//xk += deltaV;
		{
			for (int ii = 0; ii < numOVars; ii++)
				xk(ii) += deltaV(ii);

			for (int ii = 0; ii < numEConstr; ii++)
				lambdaK(ii) += deltaV(numOVars + ii);
		}

		if (deltaV.squaredNorm() < epsilon2 * epsilon2)
		{
			break;
		}
	}

	*result = xk;

}

//-----------------------------------------------------------------------------------------------------------//
void SQP1(const CD2Func& objectiveF, const EVector& x0, const CD2Func& econstrF, const LagrangeMultMethodParams& params, EVector* result)
{
	const int numOVars = x0.rows();

	// L = F(x) + lamda * C(x)
	// gradient(L(x, lamda)) = 0  <==> gradient(F(x)) + lamda*(gradient(C(x))) = 0

	// make initial lamda0 guess.
	// lambdaStar = -[gC(xstar)^t . A . gC(xstar)]^-1 . gC(xstar)^t . A . gF(xstar), where A is nonsingular matrix that is positive definite on null space of gC(xstar)^t
	// let A be I, lambda0 = -[gC(x0)^t * gC(x0)]^-1 . gC(x0)^t . gF(x0)

	EVector lambda0(1);
	{
		EVector gC0(numOVars);
		econstrF.g(x0, &gC0);
		EVector gF0(numOVars);
		objectiveF.g(x0, &gF0);

		float uu = gC0.transpose().dot(gC0);
		lambda0(0) = -(gC0.transpose().dot(gF0)) / uu;
	}

	// L func
	LScalarFunc LFunc = [&objectiveF, &econstrF](const EVector& input, const EVector& lambdas) -> float {
		return objectiveF.f(input) + econstrF.f(input) * lambdas(0);	// there's only 1 constraint.
	};

	// gradient of L func respect to x vector
	LGradientFunc gLFunc = [&objectiveF, &econstrF, numOVars](const EVector& input, const EVector& lambdas, EVector* output) {
		EVector res1(numOVars);
		EVector res2(numOVars);
		objectiveF.g(input, &res1);
		econstrF.g(input, &res2);
		*output = res1 + res2 * lambdas(0);	// there's only 1 constraint.
	};

	// hessian of L func respect to x vector
	LHessianFunc hLFunc = [&objectiveF, &econstrF, numOVars](const EVector& input, const EVector& lambdas, EMatrix* output) {
		EMatrix res1(numOVars, numOVars);
		EMatrix res2(numOVars, numOVars);
		objectiveF.h(input, &res1);
		econstrF.h(input, &res2);
		*output = res1 + res2 * lambdas(0);	// there's only 1 constraint.
	};

	L2Func LF(LFunc, gLFunc, hLFunc);
	CD1Func econstrFs(econstrF.f, econstrF.g);

	SQPFramework(LF, x0, lambda0, 1, &econstrFs, params.m_maxIter, params.m_epsilon1, params.m_epsilon2, result);
}

//-----------------------------------------------------------------------------------------------------------//
void SQP2(const CD2Func& objectiveF, const EVector& x0, const CD2Func& inconstrF, const LagrangeMultMethodParams& params, EVector* result)
{
	// change inequality constraint to equality constraint by introducing slack variable.
	// oinput is the original variables vector
	// ninput is the variables vector with slack variables.

	int numOVars = x0.rows();

	ScalarFunc nF = [numOVars, &objectiveF](const EVector& ninput) -> float {
		ASSERT(ninput.rows() == numOVars + 1);
		EVector oinput = ninput;
		oinput.conservativeResize(ninput.rows() - 1);
		return objectiveF.f(oinput);
	};

	GradientFunc gnF = [numOVars, &objectiveF](const EVector& ninput, EVector* noutput) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numORows = ninput.rows() - 1;

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EVector ooutput(numORows);
		objectiveF.g(oinput, &ooutput);
		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1);
		(*noutput)(numORows) = 0.f;
	};

	HessianFunc hnF = [numOVars, &objectiveF](const EVector& ninput, EMatrix* noutput) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numORows = ninput.rows() - 1;
		
		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EMatrix ooutput(numORows, numORows);
		objectiveF.h(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1, numORows + 1);

		for (int ii = 0; ii < numORows; ii++)
			(*noutput)(ii, numORows) = (*noutput)(numORows, ii) = 0.f;

		(*noutput)(numORows, numORows) = 0.f;
	};

	// h(x) = C(x) + slk * slk == 0
	ScalarFunc HFunc = [numOVars, &inconstrF](const EVector& ninput) -> float {
		ASSERT(ninput.rows() == numOVars + 1);
		int numORows = ninput.rows() - 1;
		float slackVar = ninput(numORows);

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);
		return inconstrF.f(oinput) + slackVar * slackVar;
	};

	// gh(x) = (gC(x), 2 * slk)
	GradientFunc gHFunc = [numOVars, &inconstrF](const EVector& ninput, EVector* noutput) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numORows = ninput.rows() - 1;
		float slackVar = ninput(numORows);

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);
		
		EVector ooutput(numORows);
		inconstrF.g(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1);
		(*noutput)(numORows) = 2.f * slackVar;
	};

	// hh(x) = ...
	HessianFunc hHFunc = [numOVars, &inconstrF](const EVector& ninput, EMatrix* noutput) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numORows = ninput.rows() - 1;

		EVector oinput = ninput;
		oinput.conservativeResize(numORows);

		EMatrix ooutput(numORows, numORows);
		inconstrF.h(oinput, &ooutput);

		*noutput = ooutput;
		noutput->conservativeResize(numORows + 1, numORows + 1);

		for (int ii = 0; ii < numORows; ii++)
			(*noutput)(ii, numORows) = (*noutput)(numORows, ii) = 0.f;

		(*noutput)(numORows, numORows) = 2.f;
	};

	// new L func with converted equality constraint
	LScalarFunc LFunc = [numOVars, &nF, &HFunc](const EVector& ninput, const EVector& lambdas) -> float {
		ASSERT(ninput.rows() == numOVars + 1);
		return nF(ninput) + HFunc(ninput) * lambdas(0); // there's only 1 constraint.
	};

	// gradient of L func respect to x vector
	LGradientFunc gLFunc = [numOVars, &gnF, &gHFunc](const EVector& ninput, const EVector& lambdas, EVector* output) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numNVars = ninput.rows();
		EVector res1(numNVars);
		EVector res2(numNVars);
		gnF(ninput, &res1);
		gHFunc(ninput, &res2);
		*output = res1 + res2 * lambdas(0);
	};

	// hessian of L func respect to x vector
	LHessianFunc hLFunc = [numOVars, &hnF, &hHFunc](const EVector& ninput, const EVector& lambdas, EMatrix* output) {
		ASSERT(ninput.rows() == numOVars + 1);
		int numNVars = ninput.rows();
		EMatrix res1(numNVars, numNVars);
		EMatrix res2(numNVars, numNVars);
		hnF(ninput, &res1);
		hHFunc(ninput, &res2);
		*output = res1 + res2 * lambdas(0);
	};

	const int numNVars = numOVars + 1;

	EVector lambda0(1);
	lambda0(0) = params.m_lamda1;

	EVector nx0 = x0;
	ChangeEVector(&nx0, numOVars + 1);
	//{
	//	float initC = _C(x0);
	//	ASSERTF(initC <= 0.f, ("initial guess x0 doesn't satisfy constraint!"));
	//	if (initC >= 0.f)
	//	{
	//		return;
	//	}
	//	xk(numOVars) = sqrtf(-initC);
	//}

	L2Func LF(LFunc, gLFunc, hLFunc);
	CD1Func econstrFs(HFunc, gHFunc);

	EVector nresult(numNVars);
	SQPFramework(LF, nx0, lambda0, 1, &econstrFs, params.m_maxIter, params.m_epsilon1, params.m_epsilon2, &nresult);

	*result = nresult;
	ChangeEVector(result, numOVars);
}


//-----------------------------------------------------------------------------------------------------------//
void SQP3(const CD2Func& objectiveF, const EVector& x0, int numInconstr, const CD2Func* inconstrFs, const LagrangeMultMethodParams& params, EVector* result)
{
	static const int kMaxNumConstrs = 32;

	ASSERT(numInconstr >= 0);
	ASSERT(numInconstr <= kMaxNumConstrs);

	numInconstr = numInconstr < kMaxNumConstrs ? numInconstr : kMaxNumConstrs;

	int numOVars = x0.rows();

	// change inequality constraints to equality constraints by introducing slack variables.
	// oinput is the original variables vector
	// ninput is the variables vector with slack variables.

	ScalarFunc nF = [numOVars, &objectiveF, numInconstr](const EVector& ninput) -> float {
		ASSERT(ninput.rows() == numOVars + numInconstr);

		EVector oinput = ninput;
		ChangeEVector(&oinput, numOVars);

		return objectiveF.f(oinput);
	};

	GradientFunc gnF = [numOVars, &objectiveF, numInconstr](const EVector& ninput, EVector* noutput) {
		ASSERT(ninput.rows() == numOVars + numInconstr);

		EVector oinput = ninput;
		ChangeEVector(&oinput, numOVars);

		EVector ooutput(numOVars);
		objectiveF.g(oinput, &ooutput);
		*noutput = ooutput;

		ChangeEVector(noutput, numOVars + numInconstr);
	};

	HessianFunc hnF = [numOVars, &objectiveF, numInconstr](const EVector& ninput, EMatrix* noutput) {
		ASSERT(ninput.rows() == numOVars + numInconstr);

		EVector oinput = ninput;
		ChangeEVector(&oinput, numOVars);

		EMatrix ooutput(numOVars, numOVars);
		objectiveF.h(oinput, &ooutput);

		*noutput = ooutput;
		ChangeEMatrix(noutput, numOVars + numInconstr);
	};

	// h(x) = C(x) + slk * slk == 0
	ScalarFunc HFuncs[kMaxNumConstrs];
	for (int iFunc = 0; iFunc < numInconstr; iFunc++)
	{
		HFuncs[iFunc] = [numOVars, iFunc, numInconstr, &inconstrFs](const EVector& ninput) -> float {
			ASSERT(ninput.rows() == numOVars + numInconstr);
			float slackVar = ninput(numOVars + iFunc);

			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			return inconstrFs[iFunc].f(oinput) + slackVar * slackVar;
		};
	}

	// gh(x) = (gC(x), 2 * slk)
	GradientFunc gHFuncs[kMaxNumConstrs]; 
	for (int iFunc = 0; iFunc < numInconstr; iFunc++)
	{
		gHFuncs[iFunc] = [numOVars, iFunc, numInconstr, &inconstrFs](const EVector& ninput, EVector* noutput) {
			ASSERT(ninput.rows() == numOVars + numInconstr);
			float slackVar = ninput(numOVars + iFunc);

			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			EVector ooutput(numOVars);
			inconstrFs[iFunc].g(oinput, &ooutput);

			*noutput = ooutput;
			ChangeEVector(noutput, numOVars + numInconstr);
			(*noutput)(numOVars + iFunc) = 2.f * slackVar;
		};
	}

	// hh(x) = ...
	HessianFunc hHFuncs[kMaxNumConstrs]; 
	for (int iFunc = 0; iFunc < numInconstr; iFunc++)
	{
		hHFuncs[iFunc] = [numOVars, iFunc, numInconstr, &inconstrFs](const EVector& ninput, EMatrix* noutput) {
			ASSERT(ninput.rows() == numOVars + numInconstr);

			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			EMatrix ooutput(numOVars, numOVars);
			inconstrFs[iFunc].h(oinput, &ooutput);

			*noutput = ooutput;
			ChangeEMatrix(noutput, numOVars + numInconstr);
			(*noutput)(numOVars + iFunc, numOVars + iFunc) = 2.f;
		};
	}

	// new L func with converted equality constraints
	LScalarFunc LFunc = [numOVars, &nF, numInconstr, &HFuncs](const EVector& ninput, const EVector& lambdas) -> float {
		ASSERT(ninput.rows() == numOVars + numInconstr);
		ASSERT(numInconstr == lambdas.rows());

		float res = nF(ninput);
		for (int ii = 0; ii < numInconstr; ii++)
			res += HFuncs[ii](ninput) * lambdas(ii);

		return res;
	};

	// gradient of L func respect to x vector
	LGradientFunc gLFunc = [numOVars, numInconstr, &gnF, &gHFuncs](const EVector& ninput, const EVector& lambdas, EVector* output) {
		ASSERT(ninput.rows() == numOVars + numInconstr);
		ASSERT(numInconstr == lambdas.rows());
		int numNVars = ninput.rows();

		EVector res1(numNVars);
		gnF(ninput, &res1);
		*output = res1;

		for (int ii = 0; ii < numInconstr; ii++)
		{
			EVector res2(numNVars);
			gHFuncs[ii](ninput, &res2);
			*output += res2 * lambdas(ii);
		}
	};

	// hessian of L func respect to x vector
	LHessianFunc hLFunc = [numOVars, numInconstr, &hnF, &hHFuncs](const EVector& ninput, const EVector& lambdas, EMatrix* output) {
		ASSERT(ninput.rows() == numOVars + numInconstr);
		ASSERT(numInconstr == lambdas.rows());
		int numNVars = ninput.rows();

		EMatrix res1(numNVars, numNVars);
		hnF(ninput, &res1);
		*output = res1;

		for (int ii = 0; ii < numInconstr; ii++)
		{
			EMatrix res2(numNVars, numNVars);
			hHFuncs[ii](ninput, &res2);

			*output += res2 * lambdas(ii);
		}
	};

	const int numNVars = numOVars + numInconstr;

	EVector lambda0(numInconstr);
	for (int ii = 0; ii < numInconstr; ii++)
		lambda0(ii) = params.m_lamda1;

	EVector nx0 = x0;
	ChangeEVector(&nx0, numNVars);
	//{
	//	float initC = _C(x0);
	//	ASSERTF(initC <= 0.f, ("initial guess x0 doesn't satisfy constraint!"));
	//	if (initC >= 0.f)
	//	{
	//		return;
	//	}
	//	xk(numOVars) = sqrtf(-initC);
	//}

	L2Func LF(LFunc, gLFunc, hLFunc);

	CD1Func econstrFs[7] = {
		CD1Func(HFuncs[0], gHFuncs[0]),
		CD1Func(HFuncs[1], gHFuncs[1]),
		CD1Func(HFuncs[2], gHFuncs[2]),
		CD1Func(HFuncs[3], gHFuncs[3]),
		CD1Func(HFuncs[4], gHFuncs[4]),
		CD1Func(HFuncs[5], gHFuncs[5]),
		CD1Func(HFuncs[6], gHFuncs[6]),
	};

	EVector nresult(numNVars);
	SQPFramework(LF, nx0, lambda0, numInconstr, econstrFs, params.m_maxIter, params.m_epsilon1, params.m_epsilon2, &nresult);

	*result = nresult;
	result->conservativeResize(numOVars);
}
