#include "../common/common_shared.h"
#include "hessian.h"
#include "scalar_vector.h"
#include "scalar_matrix.h"

//---------------------------------------------------------------------//
Hessian::Hessian()
	: m_rows(0)
	, m_matrix(nullptr)
{}

Hessian::Hessian(int n)
{
	xassert(n > 0);	
	m_rows = n;
	m_matrix = new ScalarF[n * n];
}

Hessian::~Hessian()
{
	delete[] m_matrix;
}

//---------------------------------------------------------------------//
void Hessian::Set(int row, int col, ScalarF f)
{
	xassert(row >= 0 && row < m_rows);
	xassert(col >= 0 && col < m_rows);

	m_matrix[Index(row, col)] = f;
}

void Hessian::Evaluate(const ScalarVector& input, ScalarMatrix* output) const
{
	xassert(input.GetLength() == output->GetNumRows());
	xassert(output->GetNumRows() == output->GetNumCols());

	for (int row = 0; row < m_rows; row++)
	{
		for (int col = 0; col < m_rows; col++)
		{
			float res = m_matrix[Index(row, col)](input);
			output->Set(row, col, res);
		}
	}
}

//---------------------------------------------------------------------//
//EHessian::EHessian()
//	: m_rows(0)
//	, m_matrix(nullptr)
//{}
//
//EHessian::EHessian(int n)
//{
//	xassert(n > 0);	
//	m_rows = n;
//	m_matrix = new ScalarFunc[n * n];
//}
//
//EHessian::~EHessian()
//{
//	delete[] m_matrix;
//}

//---------------------------------------------------------------------//
//void EHessian::Set(int row, int col, ScalarFunc f)
//{
//	xassert(row >= 0 && row < m_rows);
//	xassert(col >= 0 && col < m_rows);
//
//	m_matrix[Index(row, col)] = f;
//}
//
//void EHessian::Evaluate(const EVector& input, EMatrix* output) const
//{
//	xassert(input.rows() == output->rows());
//	xassert(output->rows() == output->cols());
//
//	for (int row = 0; row < m_rows; row++)
//	{
//		for (int col = 0; col < m_rows; col++)
//		{
//			float res = m_matrix[Index(row, col)](input);
//			(*output)(row, col) = res;
//		}
//	}
//}