#include "common_shared.h"
#include <Eigen/Dense>

#include "eigen_wrapper.h"
#include "bit_array.h"

//-----------------------------------------------------------------------------------------------------------//
// helper function for eigen vectors.
//-----------------------------------------------------------------------------------------------------------//
void ChangeEVector(EVector* inout, int numRows)
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

void ChangeEMatrix(EMatrix* inout, int numRows)
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

//------------------------------------------------------------------------------------------------------------//
EMatrix MatFromRowIdx(const EMatrix& A, const EVector& rowArr)
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

EMatrix MatFromColIdx(const EMatrix& A, const EVector& colArr)
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
EVector RmRows(int numRows, const EVector& rmRow)
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

EMatrix MatRmvRowIdx(const EMatrix& A, const EVector& rmRow)
{
	EVector newIndices = RmRows(A.rows(), rmRow);
	return MatFromRowIdx(A, newIndices);
}

EVector VecFromIdx(const EVector& a, const EVector& rowArr)
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

EVector VecRmvIdx(const EVector& A, const EVector& rmIdx)
{
	EVector newIndices = RmRows(A.rows(), rmIdx);
	return VecFromIdx(A, newIndices);
}

void VecReplaceRows(EVector& a, const EVector& b, const EVector& rowIdx)
{
	ASSERT(b.rows() == rowIdx.rows());
	for (int ii = 0; ii < rowIdx.rows(); ii++)
	{
		int idx = rowIdx(ii);
		a(idx) = b(ii);
	}
}

bool anyNnz(const EMatrix& A, float eps)
{
	for (int ii = 0; ii < A.rows(); ii++)
		for (int jj = 0; jj < A.cols(); jj ++)
			if (fabs(A(ii, jj)) >= eps)
				return true;
	return false;
}

// find non zeros and store them in row and col index.
void findNonzeros(const EMatrix& m, EVector* rowIdx, EVector* colIdx)
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

static bool isZero(float x, void* params)
{
	float eps = *(static_cast<float*>(params));
	return fabs(x) <= eps;
}

static bool isNnz(float x, void* params) { return !isZero(x, params); }

EVector findZeroRows(const EVector& a, float eps)
{
	return findRows(a, isZero, &eps);
}

EVector findNnzRows(const EVector& a, float eps)
{
	return findRows(a, isNnz, &eps);
}

EVector findRows(const EVector& a, condf f, void* params)
{
	static const int kMaxNumBits = 1024;
	U64 blocks[16];
	ASSERT(ExternalBitArray::DetermineNumBlocks(kMaxNumBits) == 16);
	ExternalBitArray indices(kMaxNumBits, blocks);

	int numRows = a.rows();
	for (int ii = 0; ii < numRows; ii++)
		if (f(a(ii), params))
			indices.SetBit(ii);

	int numNRows = indices.CountSetBits();
	int idx = 0;
	EVector res(numNRows);
	for (int ii = 0; ii < numRows; ii++)
		if (indices.IsBitSet(ii))
			res(idx++) = ii;
	return res;
}

//-------------------------------------------------------------------------------//
EVector VecCond(const EVector& a, condf t, void* params)
{
	int numRows = a.rows();
	EVector res(numRows);

	for (int ii = 0; ii < numRows; ii++)
		res(ii) = t(a(ii), params) ? 1 : 0;

	return res;
}

static bool iseqf(float a, void* params)
{
	float b = *(static_cast<float*>(params));
	return a == b;
}

static bool isgreaterthan(float a, void* params)
{
	float b = *(static_cast<float*>(params));
	return a > b;
}

static bool issmallerthan(float a, void* params)
{
	float b = *(static_cast<float*>(params));
	return a < b;
}

EVector VecEq(const EVector& a, float t)
{
	return VecCond(a, iseqf, &t);
}

EVector VecGt(const EVector& a, float t)
{
	return VecCond(a, isgreaterthan, &t);
}

EVector VecSt(const EVector& a, float t)
{
	return VecCond(a, issmallerthan, &t);
}

EVector VecAbs(const EVector& a)
{
	int numRows = a.rows();
	EVector res(numRows);

	for (int ii = 0; ii < numRows; ii++)
		res(ii) = a(ii) < 0.f ? -a(ii) : a(ii);

	return res;
}

EVector VecDivVec(const EVector& a, const EVector& b)
{
	ASSERT(a.rows() == b.rows());
	int numRows = a.rows();
	EVector res(numRows);
	for (int ii = 0; ii < numRows; ii++)
		res(ii) = a(ii) / b(ii);
	return res;
}

//-------------------------------------------------------------------------------//
EVector VecAnd(const EVector& a, const EVector& b)
{
	ASSERT(a.rows() == b.rows());
	int numRows = a.rows();
	EVector res(numRows);

	for (int ii = 0; ii < numRows; ii++)
		res(ii) = (a(ii) != 0.f) && (b(ii) != 0.f);

	return res;
}

EVector VecOr(const EVector& a, const EVector& b)
{
	ASSERT(a.rows() == b.rows());
	int numRows = a.rows();
	EVector res(numRows);

	for (int ii = 0; ii < numRows; ii++)
		res(ii) = (a(ii) != 0.f) || (b(ii) != 0.f);

	return res;
}

EVector VecNot(const EVector& a)
{
	int numRows = a.rows();
	EVector res(numRows);

	for (int ii = 0; ii < numRows; ii++)
		res(ii) = (a(ii) == 0.f);

	return res;
}

float VecMin(const EVector& a)
{
	if (a.rows() == 0)
		return 0.f;

	float minV = a(0);

	for (int ii = 0; ii < a.rows(); ii++)
		if (a(ii) < minV)
			minV = a(ii);
	return minV;
}

float VecMin2(const EVector& a, EVector& indices)
{
	if (a.rows() == 0)
		return 0.f;
	
	float minV = a(0);

	for (int ii = 1; ii < a.rows(); ii++)
		if (a(ii) < minV)
		{
			minV = a(ii);
		}

	int count = 0;
	for (int ii = 0; ii < a.rows(); ii++)
		if (a(ii) == minV)
			count++;

	int ind = 0;
	indices = EVector(count);
	for (int ii = 0; ii < a.rows(); ii++) 
		if (a(ii) == minV)
			indices(ind++) = ii;
	ASSERT(ind == count);

	return minV;
}

//-------------------------------------------------------------------------------//
EVector colon(int j, int k)
{
	if (j <= k)
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

	return EVector();
}

EVector VecAppend(const EVector& a, const EVector& b)
{
	int numRowsA = a.rows();
	int numRowsB = b.rows();
	int numNRows = numRowsA + numRowsB;
	EVector res(numNRows);
	res.block(0, 0, numRowsA, 1) = a;
	res.block(numRowsA, 0, numRowsB, 1) = b;
	return res;
}

EMatrix MatRowAppend(const EMatrix& a, const EMatrix& b)
{
	ASSERT(a.cols() == b.cols());
	int numNRows = a.rows() + b.rows();
	int numCols = a.cols();
	EMatrix res(numNRows, numCols);
	res.block(0, 0, a.rows(), numCols) = a;
	res.block(a.rows(), 0, b.rows(), numCols) = b;
	return res;
}

EMatrix MatColAppend(const EMatrix& a, const EMatrix& b)
{
	ASSERT(a.rows() == b.rows());
	int numRows = a.rows();
	int numNCols = a.cols() + b.cols();
	EMatrix res(numRows, numNCols);
	res.block(0, 0, numRows, a.cols()) = a;
	res.block(0, a.cols(), numRows, b.cols()) = b;
	return res;
}
//------------------------------------------------------------------------------------------------------------//
EMatrix MatAddScalar(const EMatrix& m, float s)
{
	EMatrix res = m;
	for (int ii = 0; ii < m.rows(); ii++)
		for (int jj = 0; jj < m.cols(); jj++)
			res(ii, jj) += s;
	return res;
}

EMatrix MatRand(int row, int col)
{
	EMatrix res(row, col);
	for (int ii = 0; ii < row; ii++)
		for (int jj = 0; jj < col; jj ++)
			res(ii, jj) = (float)rand() / RAND_MAX;
	return res;
}

//------------------------------------------------------------------------------------------------------------//
bool EigenLlt(const EMatrix& m, EMatrix* l)
{
	Eigen::LLT<EMatrix> lltOfM; 
	lltOfM.compute(m);
	bool isPosD = lltOfM.info() == Eigen::Success;	// whether H is positive definite matrix.
	if (isPosD)
		*l = lltOfM.matrixL();

	return isPosD;
}


void EigenValVec(const EMatrix& m, 
	void* _eigenval, 
	void* _eigenvec)
{
	Eigen::EigenSolver<EMatrix> eh(m);
	EMatrix::EigenvaluesReturnType* eigenval = (EMatrix::EigenvaluesReturnType*)_eigenval;
	Eigen::EigenSolver<EMatrix>::EigenvectorsType* eigenvec = (Eigen::EigenSolver<EMatrix>::EigenvectorsType*)_eigenvec;
	*eigenval = eh.eigenvalues();
	*eigenvec = eh.eigenvectors();
}

int FindMinEigenValIdx(void* eigenval, float* outMinEigenVal)
{
	EMatrix::EigenvaluesReturnType eigenvalP = *(static_cast<EMatrix::EigenvaluesReturnType*>(eigenval));

	// find minimun eigenvec column.
	int indminR = 0;
	int nn = eigenvalP.rows();
	float smallestVal = eigenvalP(0).real();
	for (int ii = 1; ii < nn; ii++)
	{
		if (eigenvalP(ii).real() < smallestVal)
		{
			smallestVal = eigenvalP(ii).real();
			indminR = ii;
		}
	}
	if (outMinEigenVal)
		*outMinEigenVal = smallestVal;
	return indminR;
}

float FindMinEigenVal(const EMatrix& m)
{
	Eigen::EigenSolver<EMatrix> eh(m);
	EMatrix::EigenvaluesReturnType eigenvalP = eh.eigenvalues();

	int nn = eigenvalP.rows();
	float smallEigVal = eigenvalP(0).real();
	for (int ii = 1; ii < nn; ii++)
	{
		if (eigenvalP(ii).real() < smallEigVal)
		{
			smallEigVal = eigenvalP(ii).real();
		}
	}
	return smallEigVal;
}

//-----------------------------------------------------------------------------------------------------------------------------//

void EigenQrDecomp(const EMatrix& m, EMatrix* q, EMatrix* r, EMatrix* p)
{
	Eigen::ColPivHouseholderQR<EMatrix> qrOfM(m);
	if (q != nullptr)
		*q = qrOfM.householderQ();
	if (r != nullptr)
		*r = qrOfM.matrixQR().triangularView<Eigen::Upper>();
	if (p != nullptr)
		*p = qrOfM.colsPermutation();
}

//-----------------------------------------------------------------------------------------------------------------------------//

EVector EigenColPivQrSolve(const EMatrix& A, const EVector& b)
{
	// solve A.x = b
	ASSERT(A.rows() == b.rows());
	EVector x = A.colPivHouseholderQr().solve(b);
	return x;
}

//-----------------------------------------------------------------------------------------------------------------------------//
EMatrix PseInv(const EMatrix& m, float toler)
{
	Eigen::JacobiSVD<EMatrix> svd(m, Eigen::ComputeThinU|Eigen::ComputeThinV);
	EMatrix s = svd.singularValues();
	EMatrix u = svd.matrixU();
	EMatrix v = svd.matrixV();

	ASSERT(s.cols() == 1);
	EMatrix sinv(s.rows(), s.cols());
	for (int ii = 0; ii < s.rows(); ii++)
	{
		for (int jj = 0; jj < s.cols(); jj++)
		{
			if (s(ii, jj) > toler)
				sinv(ii, jj) = 1.f / s(ii);
			else
				sinv(ii, jj) = 0.f;
		}
	}
	
	EMatrix res = v * sinv.asDiagonal() * u.transpose();
	return res;
}
