
#include "scalar_vector.h"
#include "scalar_matrix.h"

#ifndef _HESSIAN_H_
#define _HESSIAN_H_

//class EHessian
//{
//public:
//	EHessian();
//	EHessian(int n);
//	~EHessian();
//
//	void Set(int row, int col, ScalarFunc f);
//	void Evaluate(const EVector& input, EMatrix* output) const;
//
//	int Index(int r, int c) const { return r*m_rows + c; }
//
//private:
//	ScalarFunc*	m_matrix;
//	int			m_rows;		// hessian is a square matrix.
//};

typedef void (*HessianFunc)(const EVector& input, EMatrix* output);

#endif