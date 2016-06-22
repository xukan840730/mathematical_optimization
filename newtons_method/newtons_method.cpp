#include "../common/common_shared.h"
#include "../line_search/line_search_subp.h"
#include "newtons_method.h"

#define SEARCH_DIRECTION_MIN_LENGTH 0.00001f

//-------------------------------------------------------------------------------------------------------------//
void NewtonsMethod(const ScalarF F, const Gradient* g, const Hessian* H, const NewtonsMethodParams& params,
	const ScalarVector& x1, ScalarVector* result)
{
	xassert(x1.GetLength() == result->GetLength());
	
	const int numParams = x1.GetLength();

	// for newton's method, x1 needs to be sufficiently close to optima.
	ScalarVector xk = x1;
	ScalarVector gk(numParams);
	ScalarVector gkNeg(numParams);
	ScalarMatrix Gk(numParams, numParams);
	ScalarMatrix GkInv(numParams, numParams);
	ScalarMatrix L(numParams, numParams);

	ScalarVector deltaK(numParams);
	
	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk);
		if (fk < params.m_min)
		{
			result->CopyFrom(xk);
			return;
		}

		g->Evaluate(xk, &gk);
		float norm2 = gk.Norm2();
		if (norm2 < params.m_epsilon * params.m_epsilon)
		{
			result->CopyFrom(xk);
			return;
		}

		// solve G(k) * deltaK = -g(k)
		H->Evaluate(xk, &Gk);

		float v = params.m_v;

		// if Gk is not positive definite, let Gk = (Gk + vI) until Gk becomes a positive definite
		do 
		{
			bool valid = CholeskyDecomposition(Gk, &L);
			if (valid)
			{ 
				LLtInverse(&GkInv, L);
				break;
			}
			else
			{
				Gk.AddI(v);
				v *= 2.f;
			}
		} while (true);

		VectorMult(&gkNeg, gk, -1.f);
		MatrixMult(&deltaK, GkInv, gkNeg);

		// update x(k)
		xk.Add(deltaK);
	}

	result->CopyFrom(xk);
}

void QuasiNewtonSR1(const ScalarF F, const Gradient* g, const NewtonsMethodParams& params, 
	const ScalarVector& x1, ScalarVector* result)
{
	xassert(x1.GetLength() == result->GetLength());

	const int numParams = x1.GetLength();

	ScalarMatrix Hk(numParams, numParams);
	Hk.Identity();

	ScalarVector xk = x1;
	ScalarVector xk1(numParams);	// x(k+1)
	ScalarVector gk(numParams);		// g(k)
	ScalarVector gk1(numParams);	// g(k+1)
	ScalarVector sk(numParams);

	ScalarVector deltaK(numParams);
	ScalarVector gammaK(numParams);

	ScalarVector t0(numParams);
	ScalarVector t1(numParams);
	ScalarMatrix Ek(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk);
		if (fk < params.m_min)
		{
			result->CopyFrom(xk);
			return;
		}

		// evaluate g(k)
		g->Evaluate(xk, &gk);

		{
			float norm2 = gk.Norm2();
			if (norm2 < params.m_epsilon * params.m_epsilon)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		// find line search direction s(k) = -H(k) * g(k)
		MatrixMult(&sk, Hk, gk);
		sk.Multiply(-1.f);

		{
			// safe guard that search direction is very small, which means we should stop.
			float norm2 = sk.Norm2();
			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		LineSearchParams lparams;
		lparams.fMin = params.m_min;
		lparams.rho = 0.01f;
		lparams.sigma = 0.1f;
		lparams.tau1 = 9.f;
		lparams.tau2 = 0.1f;
		lparams.tau3 = 0.5f;

		// x(k+1) = x(k) + alphaK * s(k)
		float alphaK = InexactLineSearch(F, *g, sk, xk, lparams);		
		VectorMult(&deltaK, sk, alphaK);
		VectorAdd(&xk1, xk, deltaK);

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			g->Evaluate(xk1, &gk1);
			VectorSubtract(&gammaK, gk1, gk);

			// H(k+1) = H(k) + ((deltaK - H . gammaK) . (deltaK - H . gammaK)T) / ((deltaK - H . gammaK)T . gammaK)
			{
				MatrixMult(&t0, Hk, gammaK);
				VectorSubtract(&t1, deltaK, t0);
				VectorMult(&Ek, t1, t1);
				float denom = DotProd(t1, gammaK);
				Ek.DividedBy(denom);
				Hk.Add(Ek);
			}
		}

		// update x(k)
		xk.Add(deltaK);
	}

	result->CopyFrom(xk);
}

//-------------------------------------------------------------------------------------------------------------//
void QuasiNewtonDFP(const ScalarF F, const Gradient* g, const NewtonsMethodParams& params, 
	const ScalarVector& x1, ScalarVector* result)
{
	xassert(x1.GetLength() == result->GetLength());

	const int numParams = x1.GetLength();

	ScalarMatrix Hk(numParams, numParams);
	Hk.Identity();

	ScalarVector xk = x1;
	ScalarVector xk1(numParams);	// x(k+1)
	ScalarVector gk(numParams);		// g(k)
	ScalarVector gk1(numParams);	// g(k+1)
	ScalarVector sk(numParams);

	ScalarVector deltaK(numParams);
	ScalarVector gammaK(numParams);

	ScalarMatrix t0(numParams, numParams);
	ScalarMatrix t1(numParams, numParams);
	ScalarMatrix t2(numParams, numParams);
	ScalarVector t4(numParams);

	ScalarMatrix Uk(numParams, numParams);
	ScalarMatrix Vk(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk);
		if (fk < params.m_min)
		{
			result->CopyFrom(xk);
			return;
		}

		// evaluate g(k)
		g->Evaluate(xk, &gk);

		{
			float norm2 = gk.Norm2();
			if (norm2 < params.m_epsilon * params.m_epsilon)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		// find line search direction s(k) = -H(k) * g(k)
		MatrixMult(&sk, Hk, gk);
		sk.Multiply(-1.f);

		{
			// safe guard that search direction is very small, which means we should stop.
			float norm2 = sk.Norm2();
			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		LineSearchParams lparams;
		lparams.fMin = params.m_min;
		lparams.rho = 0.01f;
		lparams.sigma = 0.1f;
		lparams.tau1 = 9.f;
		lparams.tau2 = 0.1f;
		lparams.tau3 = 0.5f;

		// x(k+1) = x(k) + alphaK * s(k)
		float alphaK = InexactLineSearch(F, *g, sk, xk, lparams);
		VectorMult(&deltaK, sk, alphaK);
		VectorAdd(&xk1, xk, deltaK);

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			g->Evaluate(xk1, &gk1);
			VectorSubtract(&gammaK, gk1, gk);

			// H(k+1) = H(k) + U + V = H(k) + (deltaK . (deltaK)T) / ((deltaK)T . gammaK) - (H(k) . gammaK . (gammaK)T . H(k)) / ((gammaK)T . H . gammaK)

			// Uk
			{
				VectorMult(&Uk, deltaK, deltaK);
				float den = DotProd(deltaK, gammaK);
				Uk.DividedBy(den);
			}

			// Vk
			{
				t0.CopyFrom(Hk);
				VectorMult(&t1, gammaK, gammaK);
				MatrixMult(&t2, t0, t1);
				MatrixMult(&Vk, t2, Hk);

				MatrixMult(&t4, Hk, gammaK);
				float den = -1.f * DotProd(gammaK, t4);
				Vk.DividedBy(den);
			}

			Hk.Add(Uk);
			Hk.Add(Vk);
		}

		// update x(k)
		xk.Add(deltaK);
	}

	result->CopyFrom(xk);
}

//-------------------------------------------------------------------------------------------------------------//
void QuasiNewtonBFGS(const ScalarF F, const Gradient* g, const NewtonsMethodParams& params, 
	const ScalarVector& x1, ScalarVector* result)
{
	xassert(x1.GetLength() == result->GetLength());

	const int numParams = x1.GetLength();

	ScalarMatrix Hk(numParams, numParams);
	Hk.Identity();

	ScalarVector xk = x1;
	ScalarVector xk1(numParams);	// x(k+1)
	ScalarVector gk(numParams);		// g(k)
	ScalarVector gk1(numParams);	// g(k+1)
	ScalarVector sk(numParams);

	ScalarVector deltaK(numParams);
	ScalarVector gammaK(numParams);

	ScalarVector u0(numParams);
	ScalarMatrix u1(numParams, numParams);
	ScalarMatrix u2(numParams, numParams);

	ScalarMatrix v0(numParams, numParams);
	ScalarMatrix v1(numParams, numParams);
	ScalarMatrix v2(numParams, numParams);

	ScalarMatrix Uk(numParams, numParams);
	ScalarMatrix Vk(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk);
		if (fk < params.m_min)
		{
			result->CopyFrom(xk);
			return;
		}

		// evaluate g(k)
		g->Evaluate(xk, &gk);

		{
			float norm2 = gk.Norm2();
			if (norm2 < params.m_epsilon * params.m_epsilon)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		// find line search direction s(k) = -H(k) * g(k)
		MatrixMult(&sk, Hk, gk);
		sk.Multiply(-1.f);

		{
			// safe guard that search direction is very small, which means we should stop.
			float norm2 = sk.Norm2();
			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
			{
				result->CopyFrom(xk);
				return;
			}
		}

		LineSearchParams lparams;
		lparams.fMin = params.m_min;
		lparams.rho = 0.01f;
		lparams.sigma = 0.1f;
		lparams.tau1 = 9.f;
		lparams.tau2 = 0.1f;
		lparams.tau3 = 0.5f;

		// x(k+1) = x(k) + alphaK * s(k)
		float alphaK = InexactLineSearch(F, *g, sk, xk, lparams);
		VectorMult(&deltaK, sk, alphaK);
		VectorAdd(&xk1, xk, deltaK);

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			g->Evaluate(xk1, &gk1);
			VectorSubtract(&gammaK, gk1, gk);

			//// H(k+1) = H(k) + U + V = H(k) + (deltaK . (deltaK)T) / ((deltaK)T . gammaK) - (H(k) . gammaK . (gammaK)T . H(k)) / ((gammaK)T . H . gammaK)

			// Uk
			{
				MatrixMult(&u0, Hk, gammaK);
				float u1 = DotProd(gammaK, u0);
				float den = DotProd(deltaK, gammaK);
				float u3 = (1.f + u1 / den);

				VectorMult(&Uk, deltaK, deltaK);
				Uk.Multiply(u3);
				Uk.DividedBy(den);
			}

			// Vk
			{
				VectorMult(&v0, deltaK, gammaK);
				MatrixMult(&Vk, v0, Hk);

				VectorMult(&v1, gammaK, deltaK);
				MatrixMult(&v2, Hk, v1);

				Vk.Add(v2);

				float den = -1.f * DotProd(deltaK, gammaK);
				Vk.DividedBy(den);
			}

			Hk.Add(Uk);
			Hk.Add(Vk);
		}

		// update x(k)
		xk.Add(deltaK);
	}

	result->CopyFrom(xk);
}
