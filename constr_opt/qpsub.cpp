#include "../common/common_shared.h"
#include "../common/bit_array.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"

//static EMatrix MatrixFromRowIdx(const EMatrix& A, const ExternalBitArray& indices)
//{
//	int numNRows = indices.CountSetBits();
//	EMatrix res(numNRows, A.cols());
//
//	int rowIdx = 0;
//	for (int ii = 0; ii < indices.GetMaxBitCount(); ii++)
//	{
//		if (indices.IsBitSet(ii))
//		{
//			res.col(rowIdx) = A.row(ii);
//			rowIdx++;
//		}
//	}
//	ASSERT(rowIdx == numNRows);
//	return res;
//}

//static EMatrix MatrixFromColIdx(const EMatrix& A, const ExternalBitArray& indices)
//{
//	int numNCols = indices.CountSetBits();
//	EMatrix res(A.rows(), numNCols);
//
//	int colIdx = 0;
//	for (int ii = 0; ii < indices.GetMaxBitCount(); ii++)
//	{
//		if (indices.IsBitSet(ii))
//		{
//			res.col(colIdx) = A.col(ii);
//			colIdx++;
//		}
//	}
//	ASSERT(colIdx == numNCols);
//	return res;
//}

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

// EQNSOLV Helper function for QPSUB.
//   Finds a feasible point with respect to the equality constraints.
//   If the equalities are dependent but not consistent, warning
//   messages are given. If the equalities are dependent but consistent, 
//   the redundant constraints are removed and the corresponding variables 
//   adjusted.
void eqnsolv(const EMatrix& A, const EVector& eqix, int numVars, float eps)
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
}
