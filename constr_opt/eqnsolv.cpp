#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"
#include "eqnsolv.h"

// EQNSOLV Helper function for QPSUB.
//   Finds a feasible point with respect to the equality constraints.
//   If the equalities are dependent but not consistent, warning
//   messages are given. If the equalities are dependent but consistent, 
//   the redundant constraints are removed and the corresponding variables 
//   adjusted.
eqnres eqnsolv(EMatrix& A, EVector& b, EVector& eqix, int numVars, float eps)
{
	eqnres ret;
	ret.exitFlag = 0;

	// set tolerances
	float tolDep = 100.f * numVars * eps;

	int numEcstr = eqix.rows();
	EMatrix Aeq = MatFromRowIdx(A, eqix);
	EVector beq = VecFromIdx(b, eqix);

	// See if the equalities from a consistent system:
	//   QR factorization of A
	EMatrix Qa, Ra;
	EigenQrDecomp(Aeq, &Qa, &Ra);
	// Now need to check which is dependent
	EVector depIdx;
	{
		EVector Rdiag = Ra.diagonal();
		depIdx = findZeroRows(Rdiag, tolDep);
	}

	if (numEcstr > numVars)
	{
		depIdx = VecAppend(depIdx, colon(numVars, numEcstr - 1));
	}

	bool notConsist = false;
	if (depIdx.rows() > 0)
	{
		EMatrix t1 = MatFromColIdx(Qa, depIdx);
		EMatrix t2 = t1.transpose() * VecFromIdx(b, eqix);
		notConsist = anyNnz(t2, tolDep);

		if (notConsist)
		{
			// equality constraints are inconsistent
			ret.exitFlag = -1;
			return ret;
		}
		else
		{
			// equality constraints are consistent
			int numDepend = depIdx.nonZeros();

			// delete the redudant constraints:
			// By QR factoring the transpose, we see which columns of A'
			//   (rows of A) move to the end
			EMatrix Qat, Rat, Pat;
			EigenQrDecomp(Aeq.transpose(), &Qat, &Rat, &Pat);

			EMatrix lhs = Aeq.transpose() * Pat;
			EMatrix rhs = Qat * Rat;
			
			EVector rmIdx;
			{
				EVector i, j;
				findNonzeros(Pat, &i, &j);
				rmIdx = VecFromIdx(i, depIdx);
			}
			if (rmIdx.rows() > 0)
			{
				// remove duplicated equality constraints.
				numEcstr -= rmIdx.rows();
				A = MatRmvRowIdx(A, rmIdx);
				b = VecRmvIdx(b, rmIdx);
				eqix = colon(0, numEcstr - 1);
				ASSERT(A.rows() == b.rows());
			}
			ret.rmvIdx = rmIdx;
			ret.exitFlag = 1;
		} // consistency check
	} // dependency check

	return ret;
}
