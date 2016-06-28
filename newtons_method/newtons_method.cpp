#include "../common/common_shared.h"
#include "../line_search/line_search_subp.h"
#include "newtons_method.h"

#define SEARCH_DIRECTION_MIN_LENGTH 0.00001f

//-------------------------------------------------------------------------------------------------------------//
//void NewtonsMethod(const ScalarFunc F, const GradientFunc g, const HessianFunc H, const NewtonsMethodParams& params,
//	const EVector& x1, EVector* result)
//{
//	xassert(x1.rows() == result->rows());
//
//	const int numParams = x1.rows();
//
//	// for newton's method, x1 needs to be sufficiently close to optima.
//	EVector xk = x1;
//	EVector gk(numParams);
//	EMatrix Gk(numParams, numParams);
//	EMatrix GkInv(numParams, numParams);
//
//	EVector deltaK(numParams);
//
//	for (int iter = 1; iter <= params.m_maxIter; iter++)
//	{
//		const float fk = F(xk);
//		if (fk < params.m_min)
//		{
//			*result = xk;
//			return;
//		}
//
//		g(xk, &gk);
//		float norm2 = Norm2(gk);
//		if (norm2 < params.m_epsilon * params.m_epsilon)
//		{
//			*result = xk;
//			return;
//		}
//
//		// solve G(k) * deltaK = -g(k)
//		H(xk, &Gk);
//
//		float v = params.m_v;
//
//		// if Gk is not positive definite, let Gk = (Gk + vI) until Gk becomes a positive definite
//		do 
//		{
//			//Eigen::LLT<EMatrix> LltOfGk(Gk);
//			EMatrix L(numParams, numParams);
//
//			// TODO: replaced by Eigen library LLT.
//			bool valid = CholeskyDecomposition(Gk, &L);
//			if (valid)
//			//if (LltOfGk.info() == Eigen::NumericalIssue)
//			{
//				LLtInverse(&GkInv, L);
//				break;
//			}
//			else
//			{
//				Gk += v * EMatrix::Identity(numParams, numParams);
//				v *= 2.f;
//			}
//		} while (true);
//
//		EVector gkNeg = gk * -1.f;
//		//MatrixMult(&deltaK, GkInv, gkNeg);
//		deltaK = GkInv * gkNeg;
//
//		// update x(k)
//		xk += deltaK;
//	}
//
//	*result = xk;
//}

//void QuasiNewtonSR1(const ScalarF F, const Gradient* g, const NewtonsMethodParams& params, 
//	const ScalarVector& x1, ScalarVector* result)
//{
//	xassert(x1.GetLength() == result->GetLength());
//
//	const int numParams = x1.GetLength();
//
//	ScalarMatrix Hk(numParams, numParams);
//	Hk.Identity();
//
//	ScalarVector xk = x1;
//	ScalarVector xk1(numParams);	// x(k+1)
//	ScalarVector gk(numParams);		// g(k)
//	ScalarVector gk1(numParams);	// g(k+1)
//	ScalarVector sk(numParams);
//
//	ScalarVector deltaK(numParams);
//	ScalarVector gammaK(numParams);
//
//	ScalarVector t0(numParams);
//	ScalarVector t1(numParams);
//	ScalarMatrix Ek(numParams, numParams);
//
//	for (int iter = 1; iter <= params.m_maxIter; iter++)
//	{
//		const float fk = F(xk);
//		if (fk < params.m_min)
//		{
//			result->CopyFrom(xk);
//			return;
//		}
//
//		// evaluate g(k)
//		g->Evaluate(xk, &gk);
//
//		{
//			float norm2 = gk.Norm2();
//			if (norm2 < params.m_epsilon * params.m_epsilon)
//			{
//				result->CopyFrom(xk);
//				return;
//			}
//		}
//
//		// find line search direction s(k) = -H(k) * g(k)
//		MatrixMult(&sk, Hk, gk);
//		sk.Multiply(-1.f);
//
//		{
//			// safe guard that search direction is very small, which means we should stop.
//			float norm2 = sk.Norm2();
//			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
//			{
//				result->CopyFrom(xk);
//				return;
//			}
//		}
//
//		LineSearchParams lparams;
//		lparams.fMin = params.m_min;
//		lparams.rho = 0.01f;
//		lparams.sigma = 0.1f;
//		lparams.tau1 = 9.f;
//		lparams.tau2 = 0.1f;
//		lparams.tau3 = 0.5f;
//
//		// x(k+1) = x(k) + alphaK * s(k)
//		float alphaK = InexactLineSearch(F, *g, sk, xk, lparams);		
//		VectorMult(&deltaK, sk, alphaK);
//		VectorAdd(&xk1, xk, deltaK);
//
//		// update H(k) giving H(k+1)
//		{
//			// evaluate g(k+1)
//			g->Evaluate(xk1, &gk1);
//			VectorSubtract(&gammaK, gk1, gk);
//
//			// H(k+1) = H(k) + ((deltaK - H . gammaK) . (deltaK - H . gammaK)T) / ((deltaK - H . gammaK)T . gammaK)
//			{
//				MatrixMult(&t0, Hk, gammaK);
//				VectorSubtract(&t1, deltaK, t0);
//				VectorMult(&Ek, t1, t1);
//				float denom = DotProd(t1, gammaK);
//				Ek.DividedBy(denom);
//				Hk.Add(Ek);
//			}
//		}
//
//		// update x(k)
//		xk.Add(deltaK);
//	}
//
//	result->CopyFrom(xk);
//}

//-------------------------------------------------------------------------------------------------------------//
void QuasiNewtonDFP(const ScalarFunc F, const GradientFunc g, const NewtonsMethodParams& params, 
	const EVector& x1, void* pUserData, void* pReserved, EVector* result)
{
	xassert(x1.rows() == result->rows());

	const int numParams = x1.rows();

	EMatrix Hk(numParams, numParams);
	Hk.setIdentity();

	EVector xk = x1;
	EVector xk1(numParams);	// x(k+1)
	EVector gk(numParams);		// g(k)
	EVector gk1(numParams);	// g(k+1)

	EVector deltaK(numParams);
	EVector gammaK(numParams);

	EMatrix Uk(numParams, numParams);
	EMatrix Vk(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk, pUserData, pReserved);
		if (fk < params.m_min)
		{
			*result = xk;
			return;
		}

		// evaluate g(k)
		g(xk, &gk, pUserData, pReserved);

		{
			float norm2 = gk.squaredNorm();
			if (norm2 < params.m_epsilon * params.m_epsilon)
			{
				*result = xk;
				return;
			}
		}

		// find line search direction s(k) = -H(k) * g(k)
		EVector sk = -1.f * Hk * gk;

		{
			// safe guard that search direction is very small, which means we should stop.
			float norm2 = sk.squaredNorm();
			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
			{
				*result = xk;
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
		float alphaK = InexactLineSearch(F, *g, sk, xk, pUserData, pReserved, lparams);
		deltaK = sk * alphaK;
		xk1 = xk + deltaK;

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			g(xk1, &gk1, pUserData, pReserved);
			gammaK = gk1 - gk;

			// H(k+1) = H(k) + U + V = H(k) + (deltaK . (deltaK)T) / ((deltaK)T . gammaK) - (H(k) . gammaK . (gammaK)T . H(k)) / ((gammaK)T . H . gammaK)

			// Uk
			{
				Uk = deltaK * deltaK.transpose();
				float den = deltaK.dot(gammaK);
				Uk /= den;
			}

			// Vk
			{
				EMatrix t1 = gammaK * gammaK.transpose();
				EMatrix t2 = Hk * t1;
				Vk = t2 * Hk;

				EVector t4 = Hk * gammaK;
				float den = -1.f * gammaK.dot(t4);
				Vk /= den;
			}

			Hk += Uk;
			Hk += Vk;
		}

		// update x(k)
		xk += deltaK;
	}

	*result = xk;
}

//-------------------------------------------------------------------------------------------------------------//
void QuasiNewtonBFGS(const ScalarFunc F, const GradientFunc g, const NewtonsMethodParams& params,
	const EVector& x1, void* pUserData, void* pReserved, EVector* result)
{
	xassert(x1.rows() == result->rows());

	const int numParams = x1.rows();

	EMatrix Hk(numParams, numParams);
	Hk.setIdentity();

	EVector xk = x1;
	EVector xk1(numParams);	// x(k+1)
	EVector gk(numParams);	// g(k)
	EVector gk1(numParams);	// g(k+1)
	EVector sk(numParams);

	EMatrix Uk(numParams, numParams);
	EMatrix Vk(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = F(xk, pUserData, pReserved);
		if (fk < params.m_min)
		{
			*result = xk;
			return;
		}

		// evaluate g(k)
		g(xk, &gk, pUserData, pReserved);

		{
			float norm2 = Norm2(gk);
			if (norm2 < params.m_epsilon * params.m_epsilon)
			{
				*result = xk;
				return;
			}
		}

		// find line search direction s(k) = -H(k) * g(k)
		sk = (Hk * gk) * -1.f;

		{
			// safe guard that search direction is very small, which means we should stop.
			float norm2 = Norm2(sk);
			if (norm2 < SEARCH_DIRECTION_MIN_LENGTH * SEARCH_DIRECTION_MIN_LENGTH)
			{
				*result = xk;
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
		float alphaK = InexactLineSearch(F, *g, sk, xk, pUserData, pReserved, lparams);
		EVector deltaK = sk * alphaK;
		xk1 = xk + deltaK;

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			g(xk1, &gk1, pUserData, pReserved);
			EVector gammaK = gk1 - gk;

			//// H(k+1) = H(k) + U + V = H(k) + (deltaK . (deltaK)T) / ((deltaK)T . gammaK) - (H(k) . gammaK . (gammaK)T . H(k)) / ((gammaK)T . H . gammaK)

			// Uk
			{
				EVector u0 = Hk * gammaK;
				float u1 = gammaK.dot(u0);
				float den = deltaK.dot(gammaK);
				float u3 = (1.f + u1 / den);

				Uk = deltaK * deltaK.transpose();
				Uk *= u3 / den;
			}

			// Vk
			{
				EMatrix v0 = deltaK * gammaK.transpose();
				Vk = v0 * Hk;

				EMatrix v1 = gammaK * deltaK.transpose();
				EMatrix v2 = Hk * v1;

				Vk += v2;

				float den = -1.f * deltaK.dot(gammaK);
				Vk /= (den);
			}

			Hk += Uk;
			Hk += Vk;
		}

		// update x(k)
		xk += deltaK;
	}

	*result = xk;
}
