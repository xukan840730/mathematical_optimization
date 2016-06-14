#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "../linear_algebra/scalar_matrix.h"
#include "../line_search/line_search_subp.h"


//float func(float x1, float x2)
//{
//	float A = (x2 - x1 * x1);
//	float B = (1 - x1);
//	return 100 * A * A + B * B;
//}
//
//float fdev1(float x1, float x2)
//{
//	return -400 * (x2 - x1 * x1) * x1 - 2 * (1 - x1);
//}
//
//float fdev2(float x1, float x2)
//{
//	return 200 * (x2 - x1 * x1);
//}

float falpha(float a)
{
	return 100 * a * a * a * a + (1 - a) * (1 - a);
}

float f_alpha_dev(float a)
{
	return 400 * a * a * a - 2 * (1 - a);
}

float halpha(float a)
{
	return cosf(a) - sinf(a);
}

float h_alpha_dev(float a)
{
	return -sinf(a) - cosf(a);
}

int main()
{

	{
		FParams params;
		params.fMin = -1.f;
		params.rho = 0.01f;
		params.sigma = 0.1f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		float finalA = 0.f;

		BracketRes bracketRes = bracketing(falpha, f_alpha_dev, 0.f, 0.1f, params);

		if (bracketRes.t)
		{
			finalA = bracketRes.alpha;
			printf("done!\n");
		}
		else if (bracketRes.tB)
		{
			finalA = sectioning(falpha, f_alpha_dev, bracketRes.interval, params);
			printf("done!\n");
		}
		else
		{
			// unknown.
			//assert(false);
		}
	}

	{
		FParams params;
		params.fMin = -2.f;
		params.rho = 0.01f;
		params.sigma = 0.1f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		float finalA = 0.f;

		BracketRes bracketRes = bracketing(halpha, h_alpha_dev, 0.f, 3.f, params);

		if (bracketRes.t)
		{
			finalA = bracketRes.alpha;
			printf("done!\n");
		}
		else if (bracketRes.tB)
		{
			finalA = sectioning(halpha, h_alpha_dev, bracketRes.interval, params);
			printf("done!\n");
		}
		else
		{
			// unknown.
			//assert(false);
		}
	}

	{
		ScalarMatrix A(3, 3);
		ScalarMatrix L(3, 3);
		ScalarMatrix U(3, 3);

		ScalarMatrix tA(3, 3);
		ScalarMatrix iA(3, 3);
		{
			A.Set(0, 0, 25.f);
			A.Set(0, 1, 5.f);
			A.Set(0, 2, 1.f);

			A.Set(1, 0, 64.f);
			A.Set(1, 1, 8.f);
			A.Set(1, 2, 1.f);

			A.Set(2, 0, 144.f);
			A.Set(2, 1, 12.f);
			A.Set(2, 2, 1.f);
		}

		{
			L.Set(0, 0, 1.f);
			L.Set(0, 1, 0.f);
			L.Set(0, 2, 0.f);

			L.Set(1, 0, 2.56f);
			L.Set(1, 1, 1.f);
			L.Set(1, 2, 0.f);

			L.Set(2, 0, 5.76f);
			L.Set(2, 1, 3.5f);
			L.Set(2, 2, 1.f);
		}

		{
			U.Set(0, 0, 25.f);
			U.Set(0, 1, 5.f);
			U.Set(0, 2, 1.f);

			U.Set(1, 0, 0.f);
			U.Set(1, 1, -4.8f);
			U.Set(1, 2, -1.56f);

			U.Set(2, 0, 0.f);
			U.Set(2, 1, 0.f);
			U.Set(2, 2, 0.7f);
		}

		MatrixMult(&tA, L, U);

		LUInverse(&iA, L, U);

		ScalarMatrix ttA(3, 3);
		MatrixMult(&ttA, A, iA);

		printf("done!\n");
	}


	return 0;
}