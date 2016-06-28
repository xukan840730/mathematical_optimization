
#include "scalar_vector.h"

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

typedef void (*GradientFunc)(const EVector& input, EVector* output, void* pUserData, void* pReserved);

//class EGradient
//{
//public:
//	EGradient();
//	EGradient(int n);
//	~EGradient();
//
//	void Set(int i, ScalarFunc f);
//	void Evaluate(const EVector& input, EVector* output) const;
//
//	int GetLength() const { return m_length; }
//
//private:
//	ScalarFunc*	m_funcs;
//	int			m_length;
//};

#endif