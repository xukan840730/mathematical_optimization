#include "../common/common_shared.h"
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
