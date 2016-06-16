#include "../common/common_shared.h"
#include "scalar_matrix.h"

//----------------------------------------------------------------------------------------------------//
ScalarMatrix::ScalarMatrix(int r, int c)
{
	xassert(r > 0 && c > 0);
	m_rows = r;
	m_cols = c;
	// TODO: replace new with my own pool allocator.
	m_matrix = new float[GetSize()];
}

ScalarMatrix::ScalarMatrix(const ScalarMatrix& m)
{
	xassert(m.m_rows > 0 && m.m_cols > 0);
	m_rows = m.m_rows;
	m_cols = m.m_cols;
	m_matrix = new float[GetSize()];
	memcpy(m_matrix, m.m_matrix, GetSize() * sizeof(float));
}

ScalarMatrix::~ScalarMatrix()
{
	delete [] m_matrix;
}

void ScalarMatrix::Identity()
{
	xassert(m_rows > 0 && m_cols > 0);

	for (int row = 0; row < m_rows; row++)
	{
		for (int col = 0; col < m_cols; col++)
		{
			m_matrix[Index(row, col)] = (row == col) ? 1.f : 0.f;
		}
	}
}

void ScalarMatrix::CopyFrom(const ScalarMatrix& m)
{
	xassert(m_rows == m.m_rows);
	xassert(m_cols == m.m_cols);

	memcpy(m_matrix, m.m_matrix, GetSize() * sizeof(float));
}

//----------------------------------------------------------------------------------------------------//

float ScalarMatrix::Get(int r, int c) const
{
	xassert(r >= 0 && r < m_rows);
	xassert(c >= 0 && c < m_cols);
	return m_matrix[Index(r, c)];
}

void ScalarMatrix::Set(int r, int c, float val)
{
	xassert(r >= 0 && r < m_rows);
	xassert(c >= 0 && c < m_cols);
	m_matrix[Index(r, c)] = val;
}


//----------------------------------------------------------------------------------------------------//
// multiply
//----------------------------------------------------------------------------------------------------//

void MatrixMult(ScalarMatrix* result, const ScalarMatrix& m1, const ScalarMatrix& m2)
{
	xassert(result->GetNumRows() == m1.GetNumRows() && result->GetNumCols() == m2.GetNumCols());
	xassert(m1.GetNumCols() == m2.GetNumRows());

	int sumCount = m1.GetNumCols();
	xassert(sumCount > 0);

	int resIndex = 0;

	for (int i=0; i<result->GetNumRows(); i++)
	{
		int m1IndexStart = i*m1.GetNumCols();

		for (int j=0; j<result->GetNumCols(); j++)
		{
			int m1Index = m1IndexStart;
			int m2Index = j;

			float sum = 0.0f;
			for (int k=0; k<sumCount; k++)
			{
				sum += m1.m_matrix[m1Index]*m2.m_matrix[m2Index];

				++m1Index;
				m2Index += m2.GetNumCols();
			}

			result->m_matrix[resIndex++] = sum;
		}
	}
}

void MatrixMult(ScalarVector* result, const ScalarMatrix& m, const ScalarVector& v)
{
	xassert(result->GetLength() == v.GetLength());
	xassert(v.GetLength() == m.GetNumCols());

	for (int ii = 0; ii < v.GetLength(); ii++)
	{
		float sum = 0.f;

		for (int col = 0; col < v.GetLength(); col++)
		{
			sum += v.Get(col) * m.Get(ii, col);
		}

		result->Set(ii, sum);
	}
}

//----------------------------------------------------------------------------------------------------//
// matrix inversion
//----------------------------------------------------------------------------------------------------//
void LUInverse(ScalarMatrix* result, const ScalarMatrix& L, const ScalarMatrix& U)
{
	xassert(L.GetNumRows() == U.GetNumCols());
	xassert(L.GetNumRows() == U.GetNumCols());
	xassert(result->GetNumRows() == result->GetNumCols());
	xassert(L.GetNumRows() == U.GetNumRows());
	xassert(result->GetNumRows() == U.GetNumRows());

	//  L * U * B = I, B is the result matrix.
	// let Z = U * B

	int numRows = L.GetNumRows();
	int numCols = L.GetNumCols();

	ScalarMatrix Z(numRows, numCols);

	for (int col = 0; col < numCols; col++)
	{
		// L * Z = I, L is known, solve Z by forward substitution

		for (int row = 0; row < numRows; row++)
		{
			float c = (row == col) ? 1.f : 0.f;

			float sum = 0.f;
			for (int kk = 0; kk < row; kk++)
			{
				sum += L.Get(row, kk) * Z.Get(kk, col);
			}

			float co = L.Get(row, row);
			xassert(fabsf(co) > NDI_FLT_EPSILON);
			float zz = (c - sum) / co;

			Z.Set(row, col, zz);
		}

		// U * B = Z, U is known, solve B by back substitution

		for (int row = numRows - 1; row >= 0; row--)
		{
			float c = Z.Get(row, col);

			float sum = 0.f;
			for (int kk = row + 1; kk < numRows; kk++)
			{
				sum += U.Get(row, kk) * result->Get(kk, col);
			}

			float co = U.Get(row, row);
			xassert(fabsf(co) > NDI_FLT_EPSILON);
			float bb = (c - sum) / co;

			result->Set(row, col, bb);
		}
	}
}

//----------------------------------------------------------------------------------------------------//
// LU Decomposition
//----------------------------------------------------------------------------------------------------//
void LUDecomposition(const ScalarMatrix& A, ScalarMatrix* L, ScalarMatrix* U)
{
	xassert(A.GetNumRows() == L->GetNumRows());
	xassert(A.GetNumRows() == L->GetNumCols());
	xassert(A.GetNumRows() == A.GetNumCols());	// it is not necessary, LU decomposition doesn't need a square matrix.

	// copy A to L and U.
	L->Identity();
	U->CopyFrom(A);

	for (int col = 0; col < A.GetNumCols(); col++)
	{
		for (int row = col + 1; row < A.GetNumRows(); row++)
		{
			float aa = U->Get(col, col);
			xassert(fabsf(aa) > NDI_FLT_EPSILON);
			float bb = U->Get(row, col);
			float scale = bb / aa;

			L->Set(row, col, scale);
		
			for (int kk = col; kk < A.GetNumCols(); kk++)
			{
				//a10 - a00 * (a10 / a00), a11 - a01 * (a10/ a00)
				float ee = U->Get(row, kk) - U->Get(col, kk) * scale;	
				U->Set(row, kk, ee);
			}
		}
	}
}

//----------------------------------------------------------------------------------------------------//
// matrix inversion.
//----------------------------------------------------------------------------------------------------//
void MatrixInverse(ScalarMatrix* result, const ScalarMatrix& A)
{
	int numRows = A.GetNumRows();
	int numCols = A.GetNumCols();

	ScalarMatrix L(numRows, numCols);
	ScalarMatrix U(numRows, numCols);
	
	LUDecomposition(A, &L, &U);

	LUInverse(result, L, U);
}