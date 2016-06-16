#include "../common/common_shared.h"
#include "newtons_method.h"

void NewtonsMethod(const ScalarF F, const Gradient* g, const Hessian* H, const ScalarVector& x1, ScalarVector* result)
{
	xassert(x1.GetLength() == result->GetLength());
	
	const int numParams = x1.GetLength();

	// for newton's method, x1 needs to be sufficiently close to optima.
	ScalarVector xk = x1;
	ScalarVector gk(numParams);
	ScalarVector gkNeg(numParams);
	ScalarMatrix Gk(numParams, numParams);
	ScalarMatrix GkInv(numParams, numParams);

	ScalarVector phok(numParams);
	
	for (int iter = 1; iter <= 5; iter++)
	{
		const float fk = F(xk);

		g->Evaluate(xk, &gk);
		VectorMult(&gkNeg, gk, -1.f);

		// solve G(k) * pho = -g(k)
		H->Evaluate(xk, &Gk);

		MatrixInverse(&GkInv, Gk);

		MatrixMult(&phok, GkInv, gkNeg);

		// update x(k) to x(k+1)
		xk.Add(phok);
	}

	result->CopyFrom(xk);
}