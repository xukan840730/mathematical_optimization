#include "../common/common_shared.h"
#include "augmented_lagrangian.h"
#include "../unconstrained_optimization/newtons_method.h"

//-----------------------------------------------------------------------------------------------------------//
// helper function for eigen vectors.
//-----------------------------------------------------------------------------------------------------------//
static void ChangeEVector(EVector* inout, int numRows)
{
	int numORows = inout->rows();
	if (numORows < numRows)
	{
		inout->conservativeResize(numRows);
		(*inout).block(numORows, 0, numRows - numORows, 1).setZero();
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

		(*inout).block(numORows, 0, numRows - numORows, numRows - numORows).setZero();
		(*inout).block(0, numORows, numRows - numORows, numRows - numORows).setZero();
		(*inout).block(numORows, numORows, numRows - numORows, numRows - numORows).setZero();
	}
	else
	{
		inout->conservativeResize(numRows, numRows);
	}
}

//-----------------------------------------------------------------------------------------------------------//
static OptResult ALFramework(const CD1Func& objectiveF, const EVector& x0, const EVector& lambda0, int numEConstr, const CD1Func* econstrFs, int maxIter, float epsilon2)
{
	EVector xk = x0;
	EVector lambdaK = lambda0;

	int numVars = x0.rows();
	
	// main loop.
	int iter = 1;
	for (; iter <= maxIter; iter++)
	{
		float c = 1.f;
		{
			float minErr = NDI_FLT_MAX;
			EVector err(numEConstr);
			for (int ii = 0; ii < numEConstr; ii++)
			{
				float e = (*econstrFs[ii].f)(xk);
				err(ii) = e;
				if (abs(e) < minErr)
					minErr = abs(e);
			}
			if (minErr > NDI_FLT_EPSILON && minErr < 1.f)
				c = 1 / minErr;
		}

		// let AL function be f(x) + c*(e(x) - theta)^2 / 2, and ignoring constant term.
		// ALfunc  = f(x) + c*e(x)^2/2 - c*theta*e(x), let lambda = c*theta
		// as every iteration, lambdas should be closer and closer to lagrangian multipliers
		ScalarFunc ALFunc = [c, numEConstr, &econstrFs, objectiveF, &lambdaK](const EVector& input) {
			ASSERT(numEConstr == lambdaK.rows());
			float res = (*objectiveF.f)(input);

			for (int ii = 0; ii < numEConstr; ii++)
			{
				float e = (*econstrFs[ii].f)(input);
				res += c*e*e / 2.f + lambdaK(ii)*e;
			}

			return res;
		};

		// gALfunc = gf(x) - c*theta*ge(x) + c*e(x)*ge(x), let lambda = c*theta
		// so gALfunc = gf(x) + lambda*ge(x) + c*e(x)*ge(x)
		GradientFunc gALFunc = [c, numEConstr, &econstrFs, objectiveF, &lambdaK](const EVector& input, EVector* output) {
			ASSERT(numEConstr == lambdaK.rows());

			EVector res(input.rows());
			(*objectiveF.g)(input, &res);

			for (int ii = 0; ii < numEConstr; ii++)
			{
				float e = (*econstrFs[ii].f)(input);
				EVector ge(input.rows());
				(*econstrFs[ii].g)(input, &ge);

				res += lambdaK(ii) * ge + c * e * ge;
			}

			*output = res;
		};

		{
			float fAL = ALFunc(xk);
			float f = (*objectiveF.f)(xk);
			if ((fAL - f) * (fAL - f) < epsilon2 * epsilon2)
			{
				break;
			}
		}

		// use BFGS to optimize sub problem.
		CD1Func subP(ALFunc, gALFunc);

		NewtonsMethodParams subParams;
		subParams.m_min = -NDI_FLT_MAX;
		subParams.m_maxIter = 100;

		EVector subResult(numVars);
		QuasiNewtonBFGS(subP, xk, subParams, &subResult);

		xk = subResult;
		// also update lambdaK to approximate lagrangian multipliers
		for (int ii = 0; ii < numEConstr; ii++)
		{
			lambdaK(ii) += c * (*econstrFs[ii].f)(xk);
		}
	}

	OptResult result;
	result.xstar = xk;
	result.numIter = iter;

	return result;
}

//-----------------------------------------------------------------------------------------------------------//
OptResult ALMethod(const CD1Func& objectiveF, const EVector& x0, int numEConstr, const CD1Func* econstrFs, const LagrangeMultMethodParams& params)
{
	EVector lambda0(numEConstr);
	for (int ii = 0; ii < numEConstr; ii++)
		lambda0(ii) = 0.f;

	EVector xk = x0;

	static const int kMaxNumConstrs = 32;

	//const CD1Func* pEconstrFs[kMaxNumConstrs];
	//for (int ii = 0; ii < numEConstr; ii++)
	//	pEconstrFs[ii] = &econstrFs[ii];

	return ALFramework(objectiveF, x0, lambda0, numEConstr, econstrFs, params.m_maxIter, params.m_epsilon2);
}

//-----------------------------------------------------------------------------------------------------------//
OptResult ALMethod(const CD1Func& objectiveF, const EVector& x0, int numEConstr, const CD1Func* econstrFs, int numInconstr, const CD1Func* inconstrFs, const LagrangeMultMethodParams& params)
{
	static const int kMaxNumConstrs = 32;

	const int numOVars = x0.rows();
	const int numNVars = numOVars + numInconstr;

	// change equality constraints to take extra parameter.
	ScalarFunc EFuncs[kMaxNumConstrs];
	GradientFunc gEFuncs[kMaxNumConstrs];
	for (int iFunc = 0; iFunc < numEConstr; iFunc++)
	{
		EFuncs[iFunc] = [numOVars, iFunc, numEConstr, &econstrFs](const EVector& ninput) -> float {
			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			return (*econstrFs[iFunc].f)(oinput);
		};

		gEFuncs[iFunc] = [numOVars, numNVars, iFunc, numEConstr, &econstrFs](const EVector& ninput, EVector* noutput) {
			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			EVector ooutput(numOVars);
			(*econstrFs[iFunc].g)(oinput, &ooutput);

			*noutput = ooutput;
			ChangeEVector(noutput, numNVars);
		};
	}

	// change inequality constraints to equality constraints by introducing slack variables.
	// oinput is the original variables vector
	// ninput is the variables vector with slack variables.

	// h(x) = C(x) + slk * slk == 0
	// gh(x) = (gC(x), 2 * slk)
	ScalarFunc HFuncs[kMaxNumConstrs];
	GradientFunc gHFuncs[kMaxNumConstrs];
	for (int iFunc = 0; iFunc < numInconstr; iFunc++)
	{
		HFuncs[iFunc] = [numOVars, iFunc, numInconstr, &inconstrFs](const EVector& ninput) -> float {
			ASSERT(ninput.rows() == numOVars + numInconstr);
			float slackVar = ninput(numOVars + iFunc);

			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			return (*inconstrFs[iFunc].f)(oinput) + slackVar * slackVar;
		};

		gHFuncs[iFunc] = [numOVars, iFunc, numInconstr, &inconstrFs](const EVector& ninput, EVector* noutput) {
			ASSERT(ninput.rows() == numOVars + numInconstr);
			float slackVar = ninput(numOVars + iFunc);

			EVector oinput = ninput;
			ChangeEVector(&oinput, numOVars);

			EVector ooutput(numOVars);
			(*inconstrFs[iFunc].g)(oinput, &ooutput);

			*noutput = ooutput;
			ChangeEVector(noutput, numOVars + numInconstr);
			(*noutput)(numOVars + iFunc) = 2.f * slackVar;
		};
	}

	ScalarFunc nF = [numOVars, &objectiveF, numInconstr](const EVector& ninput) -> float {
		ASSERT(ninput.rows() == numOVars + numInconstr);

		EVector oinput = ninput;
		ChangeEVector(&oinput, numOVars);

		return (*objectiveF.f)(oinput);
	};

	GradientFunc gnF = [numOVars, &objectiveF, numInconstr](const EVector& ninput, EVector* noutput) {
		ASSERT(ninput.rows() == numOVars + numInconstr);

		EVector oinput = ninput;
		ChangeEVector(&oinput, numOVars);

		EVector ooutput(numOVars);
		(*objectiveF.g)(oinput, &ooutput);
		*noutput = ooutput;

		ChangeEVector(noutput, numOVars + numInconstr);
	};

	CD1Func nobjectiveF(nF, gnF);

	int numTotalEConstr = numEConstr + numInconstr;

	// fill final function array.
	CD1Func totalEconstrFs[kMaxNumConstrs];
	for (int ii = 0; ii < numEConstr; ii++)
	{
		totalEconstrFs[ii] = CD1Func(EFuncs[ii], gEFuncs[ii]);
	}

	for (int ii = 0; ii < numInconstr; ii++)
	{
		totalEconstrFs[numEConstr + ii] = CD1Func(HFuncs[ii], gHFuncs[ii]);
	}

	// initialize the variables.
	EVector nx0 = x0;
	ChangeEVector(&nx0, numNVars);
	for (int ii = 0; ii < numInconstr; ii++)
		nx0(numOVars + ii) = 0.1f;

	EVector lambda0(numTotalEConstr);
	for (int ii = 0; ii < numTotalEConstr; ii++)
		lambda0(ii) = 0.f;

	OptResult result = ALFramework(nobjectiveF, nx0, lambda0, numTotalEConstr, totalEconstrFs, params.m_maxIter, params.m_epsilon2);
	ChangeEVector(&result.xstar, numOVars);

	return result;
}

