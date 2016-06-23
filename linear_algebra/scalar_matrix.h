#include "scalar_vector.h"

#ifndef _SCALAR_MATRIX_H_
#define _SCALAR_MATRIX_H_

class ScalarMatrix
{
	friend void MatrixMult(ScalarMatrix* result, const ScalarMatrix& m1, const ScalarMatrix& m2);

private:
	float*		m_matrix;
	int			m_rows;
	int			m_cols;

public:
	ScalarMatrix(int r, int c);
	ScalarMatrix(const ScalarMatrix& m);
	~ScalarMatrix();

	void Identity();
	void ZeroOut();

	void CopyFrom(const ScalarMatrix& m);
	void Transpose(ScalarMatrix* T) const;

	int GetNumRows() const { return m_rows; }
	int GetNumCols() const { return m_cols; }
	int GetSize() const { return m_rows * m_cols; }

	int Index(int r, int c) const { return r*m_cols + c; }
	float Get(int r, int c) const;
	void Set(int r, int c, float val);

	ScalarVector GetRow(int r) const;
	ScalarVector GetCol(int c) const;

	void Add(const ScalarMatrix& b);
	void AddI(float v);	// Add v * I.
	void Multiply(float v);
	void DividedBy(float v);
};

void VectorMult(ScalarMatrix* result, const ScalarVector& a, const ScalarVector& b);	// [n,1] * [1,n] => [n,n] matrix
void MatrixMult(ScalarMatrix* result, const ScalarMatrix& m1, const ScalarMatrix& m2);
void MatrixMult(ScalarVector* result, const ScalarMatrix& m, const ScalarVector& v);

// use LU matrix to calculate the inverse matrix.
void LUInverse(ScalarMatrix* result, const ScalarMatrix& L, const ScalarMatrix& U);

// use LLt matrix to calculate the inverse matrix. 
void LLtInverse(ScalarMatrix* result, const ScalarMatrix& L);
void LLtInverse(EMatrix* result, const EMatrix& L);

// LU decomposition
void LUDecomposition(const ScalarMatrix& A, ScalarMatrix* L, ScalarMatrix* U);

// Cholesky decomposition, decompose a symmetric positive definite matrix into L and Lt. return false if input matrix is not positive definite or symmetric.
bool CholeskyDecomposition(const ScalarMatrix& A, ScalarMatrix* L);
bool CholeskyDecomposition(const EMatrix& A, EMatrix* L);

// matrix inversion. make sure the result matrix is allocated before calling inversion.
void MatrixInverse(ScalarMatrix* result, const ScalarMatrix& A);

#endif