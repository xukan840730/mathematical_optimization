#include "../common/common_shared.h"
#include "../common/bit_array.h"
#include "../linear_algebra/scalar_matrix.h"
#include "quadratic-programming.h"
#include <Eigen/QR>
#include <Eigen/LU>

// to find a point satisfy constraints Ax <= b, x = xstart + t*(xend - xstart), and maximum t.
FeasibilityRes LConstrFeasibility(const EVector& xstart, const EVector& xend, const EMatrix& A, const EVector& b)
{
	ASSERT(xstart.rows() == xend.rows());
	ASSERT(A.rows() == b.rows());

	FeasibilityRes res;
	res.type = FeasibilityRes::kImpossible;

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

		if (abs(lhi) < NDI_FLT_EPSILON && rhi < 0.f)
		{
			// e.x. 0.00000 * t <= -1, no matter what t is, it can't be satisfied.
			res.type = FeasibilityRes::kImpossible;
			break;
		}
		else if (lhi > 0.f)
		{
			float tt = rhi / lhi;
			if (tt < maximumT)
			{
				maximumT = tt;
				vconstrIdx = ii;
			}
		}
		else if (lhi < 0.f)
		{
			// lhi is smaller than 0, rhi is greater than 0, so that t from [0-1] must satisfy the constraint.
		}
	}

	res.t = maximumT;
	res.vconstrIdx = vconstrIdx;

	if (vconstrIdx >= 0)
		res.type = FeasibilityRes::kInfeasible;
	else
		res.type = FeasibilityRes::kFeasible;

	return res;
}

static bool QuadProgNoConstr(const EMatrix& H, const EVector& q, EVector* xstar)
{
	EVector res = H.colPivHouseholderQr().solve(-q);
	bool slnEst = (H * res).isApprox(-q, 0.00001f);
	*xstar = res;
	return slnEst;
}

//------------------------------------------------------------------------------------------------------//
EQuadProgRes EQuadProg(const EMatrix& H, const EVector& q, const EMatrix& Aeq, const EVector& beq)
{
	ASSERT(Aeq.rows() <= Aeq.cols());	// otherwise the linear system is overdetermined and could not have a solution.
	ASSERT(Aeq.rows() == beq.rows());
	ASSERT(H.rows() == H.cols());
	ASSERT(H.rows() == q.rows());
	ASSERT(Aeq.cols() == H.cols());

	EQuadProgRes result;
	result.success = false;

	EMatrix AT = Aeq.transpose();

	int mm = AT.rows();
	int nn = AT.cols();
	ASSERT(mm >= nn);

	Eigen::HouseholderQR<EMatrix> qrOfA(AT);
	EMatrix qOfAT = qrOfA.householderQ();
	EMatrix rOfAT = qrOfA.matrixQR().triangularView<Eigen::Upper>();

	EMatrix qHat = qOfAT.block(0, 0, mm, nn);
	EMatrix qN = qOfAT.block(0, nn, mm, mm - nn);

	EMatrix rHat = rOfAT.block(0, 0, nn, nn);

	// solve A * delta = B;
	EMatrix rHatT = rHat.transpose();
	EMatrix uu = rHatT.colPivHouseholderQr().solve(beq);
	bool slnEst = (rHatT * uu).isApprox(beq, 0.0001f);
	if (!slnEst)
		return result;

	EVector xHat = qHat * uu;
	if (mm == nn)
	{
		// the number of active constaints equals to the number of variables, the system has only 1 solution.
		result.success = true;
		result.xstar = xHat;
		return result;
	}

	// let P = Qn^t . H . Qn, if P is positive definite, which means H is positive definite on the null space of constraint matrix A.
	// so we can get a unique solution.
	EMatrix P = qN.transpose() * H * qN;

	Eigen::LLT<EMatrix> lltOfP;
	lltOfP.compute(P);

	if (lltOfP.info() == Eigen::Success)
	{
		EMatrix L = lltOfP.matrixL();
		// TODO: replaced by Eigen library LLT.
		EMatrix PInv(L.rows(), L.cols());
		LLtInverse(&PInv, L);

		EMatrix ww = (xHat.transpose() * H * qN + q.transpose() * qN);
		EMatrix vv = -1.f * PInv * ww.transpose();

		EVector vPart = qN * vv;
		ASSERT(vPart.rows() == xHat.rows());
		EVector xstar = xHat + vPart;

		result.success = true;
		result.xstar = xstar;
	}
	else if (lltOfP.info() == Eigen::NumericalIssue)
	{
		// Matrix P is unbound in null-space of constraint matrix A.
	}

	return result;
}

//------------------------------------------------------------------------------------------//
void QuadProg(const EMatrix& H, const EVector& q, const EMatrix& A, const EVector& b)
{
	ASSERT(H.rows() == q.rows());
	ASSERT(H.cols() == A.cols());
	ASSERT(A.rows() == b.rows());
	
	int numVars = H.rows();
	int numConstrs = A.rows();

	// TODO: this should be solved by a linear system.
	EVector x0(numVars);
	for (int ii = 0; ii < numVars; ii++)
		x0(ii) = 0.f;

	static const int kMaxNumConstrs = 1024;
	U64 activeSetBlocks[16];
	ExternalBitArray activeSet;
	ASSERT(ExternalBitArray::DetermineNumBlocks(kMaxNumConstrs) == 16);
	activeSet.Init(kMaxNumConstrs, activeSetBlocks);
	activeSet.SetBit(0);
	activeSet.SetBit(1);

	int maxIter = numConstrs * 2;

	EVector xk = x0;
	//EVector xkminus1 = xk;

	for (int iter = 0; iter < maxIter; iter++)
	{
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

		// solve active set EQP
		EVector xeqp;
		if (numActiveConstrs == 0)
		{
			bool success = QuadProgNoConstr(H, q, &xeqp);
			if (!success)
				break;
		}
		else
		{	
			EQuadProgRes eqpRes = EQuadProg(H, q, Ak, bk);
			if (!eqpRes.success)
				break;
			xeqp = eqpRes.xstar;
		}
	
		EVector xkminus1 = xk;
		FeasibilityRes fres = LConstrFeasibility(xkminus1, xeqp, A, b);
		
		if (fres.type == FeasibilityRes::kFeasible)
		{
			EVector xkminus1 = xk;
			xk = xeqp;

			// calculate lagrangian multipliers after we get xeqp.
			// gF(x) + lambdaK * gC(x) = 0
			// => Hx + q + lambdaK * A^t = 0
			// => A^t * lambdaK = -(H.x + q)
			EMatrix rhm = -1.f * (H * xk + q);
			EVector lambdaK = Ak.transpose().colPivHouseholderQr().solve(rhm);
			ASSERT(lambdaK.rows() == numActiveConstrs);

			bool allPositive = true;
			int rmActConstrIdx = -1;
			float smallestLambda = 0.f;
			for (int ii = 0; ii < numActiveConstrs; ii++)
			{
				if (lambdaK(ii) < 0)	// TODO: < or <= ??
				{
					allPositive = false;
					if (lambdaK(ii) < smallestLambda)
					{
						smallestLambda = lambdaK(ii);
						rmActConstrIdx = ii;
					}
				}
			}

			if (allPositive)
				break;

			activeSet.ClearBit(rowIdx[rmActConstrIdx]);
		}
		else if (fres.type == FeasibilityRes::kInfeasible)
		{
			EVector xkminus1 = xk;
			xk = xkminus1 + fres.t * (xeqp - xkminus1);

			activeSet.SetBit(fres.vconstrIdx);
		}
		else
		{
			ASSERT(false);
		}

	}
	
}
