#include "ndlib/math/line_search_subp.h"
#include "ndlib/math/newtons_method.h"

#define SEARCH_DIRECTION_MIN_LENGTH 0.00001f

//-------------------------------------------------------------------------------------------------------------//
void NewtonsMethod(const CD2Func& objectiveF, const EVector& x0, const NewtonsMethodParams& params,	EVector* result)
{
	ASSERT(x0.rows() == result->rows());

	const int numParams = x0.rows();

	// for newton's method, x1 needs to be sufficiently close to optima.
	EVector xk = x0;
	EVector gk(numParams);
	EMatrix Gk(numParams, numParams);
	EMatrix GkInv(numParams, numParams);

	EVector deltaK(numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = objectiveF.f(xk);
		if (fk < params.m_min)
		{
			*result = xk;
			return;
		}

		objectiveF.g(xk, &gk);
		float norm2 = gk.squaredNorm();
		if (norm2 < params.m_epsilon * params.m_epsilon)
		{
			*result = xk;
			return;
		}

		// solve G(k) * deltaK = -g(k)
		objectiveF.h(xk, &Gk);

		float v = params.m_v;

		// if Gk is not positive definite, let Gk = (Gk + vI) until Gk becomes a positive definite
		do 
		{
			Eigen::LLT<EMatrix> LltOfGk;
			LltOfGk.compute(Gk);

			if (LltOfGk.info() == Eigen::Success)
			{
				EMatrix L = LltOfGk.matrixL();
				// TODO: replaced by Eigen library LLT.
				LLtInverse(&GkInv, L);
				break;
			}
			else if (LltOfGk.info() == Eigen::NumericalIssue)
			{
				Gk += v * EMatrix::Identity(numParams, numParams);
				v *= 2.f;
			}
			else
			{
				ALWAYS_ASSERTF(false, ("Never saw this before!"));
			}
		} while (true);

		EVector gkNeg = gk * -1.f;
		deltaK = GkInv * gkNeg;

		// update x(k)
		xk += deltaK;
	}

	*result = xk;
}

//-------------------------------------------------------------------------------------------------------------//
void QuasiNewtonBFGS(const CD1Func& objectiveF, const EVector& x0, const NewtonsMethodParams& params, EVector* result)
{	
	ASSERT(x0.rows() == result->rows());

	const int numParams = x0.rows();

	EMatrix Hk(numParams, numParams);
	Hk.setIdentity();

	EVector xk = x0;
	EVector xk1(numParams);	// x(k+1)
	EVector gk(numParams);	// g(k)
	EVector gk1(numParams);	// g(k+1)
	EVector sk(numParams);

	EMatrix Uk(numParams, numParams);
	EMatrix Vk(numParams, numParams);

	for (int iter = 1; iter <= params.m_maxIter; iter++)
	{
		const float fk = objectiveF.f(xk);
		if (fk < params.m_min)
		{
			*result = xk;
			return;
		}

		// evaluate g(k)
		objectiveF.g(xk, &gk);

		{
			float norm2 = gk.squaredNorm();
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
		float alphaK = InexactLineSearch(objectiveF.f, objectiveF.g, sk, xk, lparams);
		EVector deltaK = sk * alphaK;
		xk1 = xk + deltaK;

		// update H(k) giving H(k+1)
		{
			// evaluate g(k+1)
			objectiveF.g(xk1, &gk1);
			EVector gammaK = gk1 - gk;

			// H(k+1) = H(k) + U + V
			// U = (1 + (gammaKT . H . gammaK) / (deltaKT . gammaK)) * ((deltaK . deltaKT) / deltaKT . gammaK)
			// V = -1 * (deltaK . gammakT . H + H . gammaK . deltaKT) / (deltaKT . gammaK)

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
