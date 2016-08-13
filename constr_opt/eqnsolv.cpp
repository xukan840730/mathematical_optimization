#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"


// EQNSOLV Helper function for QPSUB.
//   Finds a feasible point with respect to the equality constraints.
//   If the equalities are dependent but not consistent, warning
//   messages are given. If the equalities are dependent but consistent, 
//   the redundant constraints are removed and the corresponding variables 
//   adjusted.
int eqnsolv(EMatrix& A, EVector& b, EVector& eqix, int numVars, float eps)
{
	int ret = 0;

	// set tolerances
	float tolDep = 100.f * numVars * eps;
	float tolCons = 1e-10;

	int numEcstr = eqix.rows();
	EMatrix Aeq = MatrixFromRowIdx(A, eqix);
	EVector beq = VectorFromIdx(b, eqix);

	// See if the equalities from a consistent system:
	//   QR factorization of A
	EMatrix Qa, Ra;
	EigenQrDecomp(Aeq, &Qa, &Ra);
	// Now need to check which is dependent
	EVector depIdx;
	{
		EVector Rdiag = Ra.diagonal();
		depIdx = findZeroIdx(Rdiag, tolDep);
	}

	if (numEcstr > numVars)
	{
		depIdx = VecAppend(depIdx, colon(numVars, numEcstr - 1));
	}

	bool notConsist = false;
	if (depIdx.rows() > 0)
	{
		EMatrix t1 = MatrixFromColIdx(Qa, depIdx);
		EMatrix t2 = t1.transpose() * VectorFromIdx(b, eqix);
		notConsist = anyNonzero(t2, tolDep);

		if (notConsist)
		{
			// equality constraints are inconsistent
			return -1;
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
				rmIdx = VectorFromIdx(i, depIdx);
			}
			if (rmIdx.rows() > 0)
			{
				numEcstr -= rmIdx.rows();
				A = MatrixRmvRowIdx(Aeq, rmIdx);
				b = VectorRmvIdx(beq, rmIdx);
				eqix = colon(0, numEcstr - 1);
				ASSERT(A.rows() == b.rows());
				ASSERT(A.rows() == eqix.rows());
			}
		} // consistency check
	} // dependency check

	return ret;
}
