#include "../common/common_shared.h"
#include "../common/bit_array.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"

static EMatrix MatrixFromRowIdx(const EMatrix& A, const ExternalBitArray& indices)
{
	int numNRows = indices.CountSetBits();
	EMatrix res(numNRows, A.cols());

	int rowIdx = 0;
	for (int ii = 0; ii < indices.GetMaxBitCount(); ii++)
	{
		if (indices.IsBitSet(ii))
		{
			res.col(rowIdx) = A.row(ii);
			rowIdx++;
		}
	}
	ASSERT(rowIdx == numNRows);
	return res;
}

static EMatrix MatrixFromColIdx(const EMatrix& A, const ExternalBitArray& indices)
{
	int numNCols = indices.CountSetBits();
	EMatrix res(A.rows(), numNCols);

	int colIdx = 0;
	for (int ii = 0; ii < indices.GetMaxBitCount(); ii++)
	{
		if (indices.IsBitSet(ii))
		{
			res.col(colIdx) = A.col(ii);
			colIdx++;
		}
	}
	ASSERT(colIdx == numNCols);
	return res;
}

void eqnsolv(const EMatrix& A, const ExternalBitArray& eqix, int numVars, float eps)
{
	float tolDep = 100.f * numVars * eps;
	float tolCons = 1e-10;

	EMatrix Aeq = MatrixFromRowIdx(A, eqix);

	EMatrix Qa, Ra;
	EigenQrDecomp(Aeq, &Qa, &Ra);
}
