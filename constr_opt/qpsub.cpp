#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"

static EMatrix MatrixFromRowIdx(const EMatrix& A, const EVector& rowArr)
{
	int numNRows = rowArr.rows();
	int numNCols = A.cols();
	int numORows = A.rows();

	EMatrix res(numNRows, numNCols);

	for (int ii = 0; ii < numNRows; ii++)
	{
		int idx = rowArr(ii);
		ASSERT(idx < numORows);
		res.row(ii) = A.row(idx);
	}

	return res;
}

static EMatrix MatrixFromColIdx(const EMatrix& A, const EVector& colArr)
{
	int numNCols = colArr.rows();
	int numNRows = A.rows();
	int numOCols = A.cols();

	EMatrix res(numNRows, numNCols);

	for (int ii = 0; ii < numNCols; ii++)
	{
		int idx = colArr(ii);
		ASSERT(idx < numOCols);
		res.col(ii) = A.col(idx);
	}

	return res;
}

static EVector VectorFromIdx(const EVector& a, const EVector& rowArr)
{
	int numNRows = rowArr.rows();
	int numORows = a.rows();

	EVector res(numNRows);

	for (int ii = 0; ii < numNRows; ii++)
	{
		int idx = rowArr(ii);
		ASSERT(idx < numORows);
		res(ii) = a(idx);
	}

	return res;
}

static bool anyNonzero(const EMatrix& A, float eps)
{
	for (int ii = 0; ii < A.rows(); ii++)
		for (int jj = 0; jj < A.cols(); jj ++)
			if (fabs(A(ii, jj)) >= eps)
				return true;
	return false;
}

static EVector findZeroIdx(const EVector& a, float eps)
{
	// figure out how many zeros first.
	int numZeros = 0;
	int numRows = a.rows();
	for (int ii = 0; ii < numRows; ii++)
		if (fabs(a(ii)) < eps) numZeros++;

	EVector res(numZeros);
	int idx = 0;
	for (int ii = 0; ii < numRows; ii++)
		if (fabs(a(ii)) < eps) res(idx++) = ii;

	ASSERT(idx == numZeros);
	return res;
}

static EVector findNnzIdx(const EVector& a, float eps)
{
	// figure out how many nonzero first
	int nnz = 0;
	int numRows = a.rows();
	for (int ii = 0; ii < numRows; ii++)
		if (fabs(a(ii)) >= eps) nnz++;

	EVector res(nnz);
	int idx = 0;
	for (int ii = 0; ii < numRows; ii++)
		if (fabs(a(ii)) >= eps) res(idx++) = ii;

	return res;
}

static EVector colon(int j, int k)
{
	int numRows = abs(k - j) + 1;
	EVector res(numRows);
	int inc = k >= j ? 1 : -1;
	int cur = j;
	for (int ii = 0; ii < numRows; ii++)
	{
		res(ii) = cur;
		cur += inc;
	}
	return res;
}

static EVector vappend(const EVector& a, const EVector& b)
{
	int numRowsA = a.rows();
	int numRowsB = b.rows();
	int numNRows = numRowsA + numRowsB;
	EVector res(numNRows);
	res.block(0, 0, numRowsA, 1) = a;
	res.block(numRowsA, 0, numRowsB, 1) = b;
	return res;
}

// EQNSOLV Helper function for QPSUB.
//   Finds a feasible point with respect to the equality constraints.
//   If the equalities are dependent but not consistent, warning
//   messages are given. If the equalities are dependent but consistent, 
//   the redundant constraints are removed and the corresponding variables 
//   adjusted.
void eqnsolv(const EMatrix& A, const EVector& b, const EVector& eqix, int numVars, int& neqcstr, float eps)
{
	// set tolerances
	float tolDep = 100.f * numVars * eps;
	float tolCons = 1e-10;

	EMatrix Aeq = MatrixFromRowIdx(A, eqix);

	// See if the equalities from a consistent system:
	//   QR factorization of A
	EMatrix Qa, Ra;
	EigenQrDecomp(Aeq, &Qa, &Ra);
	// Now need to check which is dependent
	EVector depInd;
	{
		EVector Rdiag = Ra.diagonal();
		depInd = findZeroIdx(Rdiag, tolDep);
	}

	if (neqcstr > numVars)
	{
		depInd = vappend(depInd, colon(numVars + 1, neqcstr));
	}

	bool notConsist = false;
	if (depInd.rows() > 0)
	{
		EMatrix t1 = MatrixFromColIdx(Qa, depInd);
		EMatrix t2 = t1.transpose() * VectorFromIdx(b, eqix);
		notConsist = anyNonzero(t2, tolDep);
	}
}
