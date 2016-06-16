
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
	
private:
	ScalarF*	m_funcs;
	int			m_length;
};

#endif