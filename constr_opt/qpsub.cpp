#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../common/bit_array.h"
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

// get row indices array after removing rows 
static EVector RmRows(int numRows, const EVector& rmRow)
{
	static const int kMaxNumBits = 1024;
	U64 blocks[16];
	ASSERT(ExternalBitArray::DetermineNumBlocks(kMaxNumBits) == 16);

	ASSERT(numRows <= kMaxNumBits);
	ExternalBitArray indices(numRows, blocks, true);

	// remove those rows
	for (int ii = 0; ii < rmRow.rows(); ii++)
		indices.ClearBit(rmRow(ii));

	int numNRows = indices.CountSetBits();
	
	// remaining rows
	EVector newRows(numNRows);
	int idx = 0;
	for (int ii = 0; ii < numRows; ii++)
		if (indices.IsBitSet(ii))
			newRows(idx++) = ii;
	ASSERT(idx == numNRows);
	return newRows;
}

static EMatrix MatrixRmvRowIdx(const EMatrix& A, const EVector& rmRow)
{
	EVector newIndices = RmRows(A.rows(), rmRow);
	return MatrixFromRowIdx(A, newIndices);
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

static EVector VectorRmvIdx(const EVector& A, const EVector& rmIdx)
{
	EVector newIndices = RmRows(A.rows(), rmIdx);
	return VectorFromIdx(A, newIndices);
}

static bool anyNonzero(const EMatrix& A, float eps)
{
	for (int ii = 0; ii < A.rows(); ii++)
		for (int jj = 0; jj < A.cols(); jj ++)
			if (fabs(A(ii, jj)) >= eps)
				return true;
	return false;
}

// find non zeros and store them in row and col index.
static void findNonzeros(const EMatrix& m, EVector* rowIdx, EVector* colIdx)
{
	int numNnz = 0;
	for (int ii = 0; ii < m.rows(); ii++)
		for (int jj = 0; jj < m.cols(); jj++)
			if (m(ii, jj) != 0)
				numNnz++;

	*rowIdx = EVector(numNnz);
	*colIdx = EVector(numNnz);
	int idx = 0;
	for (int col = 0; col < m.cols(); col++)
		for (int row = 0; row < m.rows(); row++)
			if (m(row, col) != 0)
			{
				(*rowIdx)(idx) = row;
				(*colIdx)(idx) = col;
				idx++;
			}
	ASSERT(idx == numNnz);
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
		depIdx = vappend(depIdx, colon(numVars, numEcstr - 1));
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
