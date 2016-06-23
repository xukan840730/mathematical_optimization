
#include "scalar_vector.h"

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

class Gradient
{
public:
	Gradient();
	Gradient(int n);
	~Gradient();

	void Set(int i, ScalarF f);
	void Evaluate(const ScalarVector& input, ScalarVector* output) const;

	int GetLength() const { return m_length; }
	
private:
	ScalarF*	m_funcs;
	int			m_length;
};

class EGradient
{
public:
	EGradient();
	EGradient(int n);
	~EGradient();

	void Set(int i, ScalarFunc f);
	void Evaluate(const EVector& input, EVector* output) const;

	int GetLength() const { return m_length; }

private:
	ScalarFunc*	m_funcs;
	int			m_length;
};

#endif