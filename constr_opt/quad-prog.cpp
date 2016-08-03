#include "../common/common_shared.h"
#include "../common/bit_array.h"
#include "../common/lin-equation.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"
#include "simplex.h"
#include "quad-prog.h"

static void PrintEVector(const EVector& v)
{
	printf("[");
	for (int i = 0; i < v.rows(); i++)
		if (i != v.rows() - 1)
			printf("%0.3f, ", v(i));
		else
			printf("%0.3f", v(i));
	printf("]");
}

static void PrintBitArray(const ExternalBitArray& b)
{
	printf("[");
	for (int i = 0; i < b.GetMaxBitCount(); i++)
		if (b.IsBitSet(i))
			printf("%d, ", i);
	printf("]");
}


struct FeasibilityRes
{
public:
	enum Type {
		kFeasible,		// the minimal of EQP is feasible in original QP problem.
		kInfeasible,	// the minimal of EQP violates constraints of original QP problem
		kImpossible,	// it's impossible to find t between [0-1] to move along x0 to xeqp
	};

	Type type;
	float t;
	int vconstrIdx;		// the first violated constrained index.
};

// to find a point satisfy constraints Ax <= b, x = xstart + t*(xend - xstart), and maximum t.
FeasibilityRes LConstrFeasibility(const EVector& xstart, const EVector& xend, const EMatrix& A, const EVector& b, const ExternalBitArray& activeSet)
{
	ASSERT(xstart.rows() == xend.rows());
	ASSERT(A.rows() == b.rows());

	FeasibilityRes res;
	res.type = FeasibilityRes::kFeasible;

	// new x = xstart + t * (xend - xstart), 
	// so original constraints Ax <= b can be written as A(xstart + t*(xend - xstart)) <= b
	// equals A * (xend - xstart) * t <= b - Ax0
	// since xstart is feasible, so that Ax0 <= b, so that c = b - Ax0 must be >= 0.

	const EVector xdelta = xend - xstart;
	const EVector rh = b - A * xstart;
	const EVector lh = A * xdelta;

	float maximumT = 1.f;
	int vconstrIdx = -1;

	for (int ii = 0; ii < A.rows(); ii++)
	{
		float lhi = lh(ii);
		float rhi = rh(ii);
		ASSERT(rhi >= 0.f);

		// we only care about inactive inequality constraint.
		if (activeSet.IsBitSet(ii))
			continue;

		LinEquation::Result lres = LinEquation::Solve(lhi, rhi, LinEquation::kLTE);
		if (lres.returnType == LinEquation::kNever)
		{
			// e.x. 0.00000 * t <= -1, no matter what t is, it can't be satisfied.
			ASSERT(false);
			res.type = FeasibilityRes::kImpossible;
			break;
		}
		else if (lres.returnType == LinEquation::kExists && lres.op == LinEquation::kLTE)
		{
			if (lres.number < maximumT)
			{
				maximumT = lres.number;
				vconstrIdx = ii;
			}
		}
	}

	res.t = maximumT;
	res.vconstrIdx = vconstrIdx;

	if (res.type == FeasibilityRes::kImpossible)
		return res;
	else if (vconstrIdx >= 0)
		res.type = FeasibilityRes::kInfeasible;
	else
		res.type = FeasibilityRes::kFeasible;

	return res;
}

int IsConstrAlwaysSatisfied(const EVector& x0, const EVector& p, const EVector& ai, const float bi, float* t)
{
	// test x = x0 + t * p, t from (0, +inf), does xp always satisfy constraint ai.x <= b.
	// Ai.(x0 + t*p) <= bi
	// t(Ai.p) <= bi - Ai.x0
	ASSERT(x0.rows() == p.rows());
	ASSERT(x0.rows() == ai.rows());
	
	float lhi = ai.dot(p);
	float rhi = bi - ai.dot(x0);

	LinEquation::Result lres = LinEquation::Solve(lhi, rhi, LinEquation::kLTE);
	if (lres.returnType == LinEquation::kAlways)
		return 1;	// always satisfy constraint
	else if (lres.returnType == LinEquation::kExists)
	{
		if (lres.op == LinEquation::kLTE)
		{
			if (lres.number <= 0)	// never satisfy constraint
				return -1;
			else 
			{
				*t = lres.number;
				return 0;
			}
		}
		else if (lres.op == LinEquation::kGTE)
		{
			if (lres.number <= 0)	// always satisfy constraint
				return 1;
			else
			{
				ASSERT(false);	// it shouldn't happen
				*t = lres.number;
				return 0;
			}
		}
		else
		{
			ASSERT(false);
		}
	}
	else
	{
		return -1;	// never satisfy constraint
	}
}

int IsConstrAlwaysSatisfied(const EVector& x0, const EVector& p, const EMatrix& A, const EVector& b, float* t)
{
	int res = 1; // always true
	float maximumT = NDI_FLT_MAX;

	for (int ii = 0; ii < A.rows(); ii++)
	{
		float t;
		EVector Ai = A.row(ii);
		float bi = b(ii);
		int cond = IsConstrAlwaysSatisfied(x0, p, Ai, bi, &t);
		switch (cond)
		{
		case 1: break;
		case 0: 
			res = 0;
			if (t < maximumT)
				maximumT = t;
		break;
		case -1: 
			res = -1;
		break;
		}

		if (res == -1)	// one of the constraints can be never be satisfied, quit!
			break;
	}
	*t = maximumT;
	return res;
}

static EVector SolveLambda(const EMatrix& H, const EVector& q, const EVector xk, const EMatrix& Ak)
{
	// gF(x) + lambdaK * gC(x) = 0
	// => Hx + q + lambdaK * A^t = 0
	// => A^t * lambdaK = -(H.x + q)
	ASSERT(H.rows() == q.rows());
	ASSERT(q.rows() == xk.rows());
	ASSERT(Ak.rows() > 0);
	EMatrix rhm = -1.f * (H * xk + q);
	//EVector lambdaK = Ak.transpose().colPivHouseholderQr().solve(rhm);
	EVector lambdaK = EigenColPivQrSolve(Ak.transpose(), rhm);
	//ASSERT(lambdaK.rows() == numActiveConstrs);

	return lambdaK;
}

static int FindNegativeLambda(const EVector& lambdaK)
{	
	int negLambdaIdx = -1;
	float smallestLambda = 0.f;
	for (int ii = 0; ii < lambdaK.rows(); ii++)
	{
		if (lambdaK(ii) < 0)
		{
			if (lambdaK(ii) < smallestLambda)
			{
				smallestLambda = lambdaK(ii);
				negLambdaIdx = ii;
			}
		}
	}
	return negLambdaIdx;
}

static void CalcActiveSet(const EMatrix& A, const EVector& b, const EVector& x0, ExternalBitArray& activeSet)
{
	ASSERT(A.rows() == b.rows());
	ASSERT(b.rows() <= activeSet.GetMaxBitCount());

	int numConstrs = b.rows();
	EVector c = A * x0 - b;
	for (int ii = 0; ii < numConstrs; ii++)
	{
		if (fabs(c(ii)) < NDI_FLT_EPSILON)
			activeSet.SetBit(ii);
		else
			activeSet.ClearBit(ii);
	}
}

//------------------------------------------------------------------------------------------//
int QuadProg(const EMatrix& H, const EVector& q, const EMatrix& A, const EVector& b)
{
	ASSERT(H.rows() == q.rows());
	ASSERT(H.cols() == A.cols());
	ASSERT(A.rows() == b.rows());

	int res = 0;
	
	int numVars = H.rows();
	int numConstrs = A.rows();

	EVector x0(numVars);
	int initFound = SolveInitFeasible(A, b, &x0);
	if (initFound != 0)
	{
		printf("QP: initial feasible not found!\n");
		return -1;
	}

	static const int kMaxNumConstrs = 1024;
	U64 activeSetBlocks[16];
	ExternalBitArray activeSet;
	ASSERT(ExternalBitArray::DetermineNumBlocks(kMaxNumConstrs) == 16);
	activeSet.Init(kMaxNumConstrs, activeSetBlocks);
	CalcActiveSet(A, b, x0, activeSet);

	int maxIter = numConstrs * 10;

	EVector xk = x0;
	//EVector xkminus1 = xk;
	
	const float epsilon = NDI_FLT_EPSILON;

	EMatrix HInv(H.rows(), H.cols()); // HInv is only calculated if H is positive definite. use it carefully. 
	EVector hstep(numVars); // hstep is valid only if H is not positive definite.

	bool isHPosD = false;
	{
		EMatrix lOfH;
		isHPosD = EigenLlt(H, &lOfH);
		if (isHPosD)
		{
			LLtInverse(&HInv, lOfH);
		}
		else
		{
			EMatrix::EigenvaluesReturnType eigenvalH;
			Eigen::EigenSolver<EMatrix>::EigenvectorsType eigenvecH;
			EigenValVec(H, &eigenvalH, &eigenvecH);

			// find minimun eigenvec column.
			int indminR = 0;
			{
				float smallestVal = eigenvalH(0).real();
				for (int ii = 1; ii < eigenvalH.rows(); ii++)
				{
					if (eigenvalH(ii).real() < smallestVal)
					{
						smallestVal = eigenvalH(ii).real();
						indminR = ii;
					}
				}
			}
			hstep = eigenvecH.col(indminR).real();
		}
	}

	for (int iter = 0; iter < maxIter; iter++)
	{
		printf("iter: %d, ", iter);
		// calculate gradient of each step
		EVector gk = H * xk + q;

		printf("xk: "); PrintEVector(xk); printf(", ");
		printf("gk: "); PrintEVector(gk); printf(", ");

		int numActiveConstrs = activeSet.CountSetBits();
		EMatrix Ak(numActiveConstrs, numVars);
		EVector bk(numActiveConstrs);
		static const int kMaxNumVars = 32;
		int rowIdx[kMaxNumVars];
		ASSERT(numActiveConstrs <= kMaxNumVars);
		{
			int activeIdx = 0;
			// fill Ak by active constraints
			for (int ii = 0; ii < numConstrs; ii++)
			{
				if (activeSet.IsBitSet(ii))
				{
					Ak.row(activeIdx) = A.row(ii);
					bk(activeIdx) = b(ii);
					rowIdx[activeIdx] = ii;
					activeIdx++;
				}
			}
		}

		printf("active constrs: "); PrintBitArray(activeSet); printf(", ");

		// solve active set EQP
		EVector xeqp;
		if (numActiveConstrs == 0)
		{
			if (isHPosD)
			{
				// H.x + q = 0;
				xeqp = HInv * -q;
				res = 0;
			}
			else
			{
				EVector stepp = hstep;
				ASSERT(stepp.rows() == H.cols());
				if (stepp.dot(gk) > epsilon) // make sure stepp is descent direction
					stepp *= -1.f;

				float t;
				int cond = IsConstrAlwaysSatisfied(xk, stepp, A, b, &t);
				if (cond == 1)	// constrs are always satisfied and quad func is unbounded.
				{
					res = 2;
					break;
				}
				else if (cond == -1)
				{
					ASSERT(false);
					break;
				}
				else
				{
					xeqp = xk + stepp * (t + 0.01);
					res = 1;
				}
			}
		}
		else
		{	
			EQuadProgRes eqpRes = EQuadProg(H, q, Ak, bk);
			if (eqpRes.type == EQuadProgRes::kUnbounded)
			{
				EVector stepp = eqpRes.stepp;
				if (stepp.dot(gk) > epsilon) // make sure stepp is descent direction
					stepp *= -1.f;
				xeqp = xk + stepp;
			}
			else
			{
				xeqp = eqpRes.xstar;
			}			
		}

		printf("xeqp: "); PrintEVector(xeqp); printf(", ");

		EVector xkminus1 = xk;
		EVector stepk = xeqp - xkminus1;

		bool isStepSmall = stepk.squaredNorm() < epsilon * epsilon; 
		if (isStepSmall)
		{
			// the step is tiny.
			if (numActiveConstrs == 0)
				break;

			EVector lambdaK = SolveLambda(H, q, xk, Ak);
			ASSERT(lambdaK.rows() == numActiveConstrs);

			printf("lambdaK: "); PrintEVector(lambdaK); printf(", ");

			int negLambdaIdx = FindNegativeLambda(lambdaK);
			if (negLambdaIdx < 0)
				break;

			ASSERT(negLambdaIdx < numActiveConstrs);
			activeSet.ClearBit(rowIdx[negLambdaIdx]);
		}
		else
		{
			FeasibilityRes fres = LConstrFeasibility(xkminus1, xeqp, A, b, activeSet);
			if (fres.type == FeasibilityRes::kFeasible)
			{
				xk = xeqp;

				if (numActiveConstrs > 0)
				{
					// calculate lagrangian multipliers after we get xeqp.
					EVector lambdaK = SolveLambda(H, q, xk, Ak);
					ASSERT(lambdaK.rows() == numActiveConstrs);
					
					printf("lambdaK: "); PrintEVector(lambdaK); printf(", ");

					int negLambdaIdx = FindNegativeLambda(lambdaK);
					if (negLambdaIdx < 0)
						break;

					ASSERT(negLambdaIdx < numActiveConstrs);
					activeSet.ClearBit(rowIdx[negLambdaIdx]);
				}
			}
			else if (fres.type == FeasibilityRes::kInfeasible)
			{
				xk = xkminus1 + fres.t * (xeqp - xkminus1);
				activeSet.SetBit(fres.vconstrIdx);
			}
			else
			{
				ASSERT(false);
			}
		}
		printf("\n");
	}

	printf("\nfinal: %d, ", res); PrintEVector(xk); printf("\n");
}
