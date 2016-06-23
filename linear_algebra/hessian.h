
#include "scalar_vector.h"
#include "scalar_matrix.h"

#ifndef _HESSIAN_H_
#define _HESSIAN_H_

class Hessian
{
public:
	Hessian();
	Hessian(int n);
	~Hessian();

	void Set(int row, int col, ScalarF f);
	void Evaluate(const ScalarVector& input, ScalarMatrix* output) const;

	int Index(int r, int c) const { return r*m_rows + c; }

private:
	ScalarF*	m_matrix;
	int			m_rows;		// hessian is a square matrix.
};

class EHessian
{
public:
	EHessian();
	EHessian(int n);
	~EHessian();

	void Set(int row, int col, ScalarFunc f);
	void Evaluate(const EVector& input, EMatrix* output) const;

	int Index(int r, int c) const { return r*m_rows + c; }

private:
	ScalarFunc*	m_matrix;
	int			m_rows;		// hessian is a square matrix.
};

#endif