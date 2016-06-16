#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "../common/common_shared.h"
#include "../linear_algebra/scalar_matrix.h"
#include "../line_search/line_search_subp.h"
#include "../newtons_method/newtons_method.h"


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

float f_test_newton(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return x1 * x1 * x1 * x1 + x1 * x2 + (1 + x2) * (1 + x2);
}

float f_test_newton_g1(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return 4 * x1 * x1 * x1 + x2;
}

float f_test_newton_g2(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return x1 + 2 * (1 + x2);
}

float f_test_newton_H00(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return 12 * x1 * x1;
}

float f_test_newton_H01(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return 1;
}

float f_test_newton_H10(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return 1;
}

float f_test_newton_H11(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	return 2;
}

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
		LineSearchParams params;
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
		LineSearchParams params;
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
		ScalarMatrix iA(3, 3);

		ScalarMatrix ttA(3, 3);

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

		// LU decomposition
		{
			ScalarMatrix L(3, 3);
			ScalarMatrix U(3, 3);

			LUDecomposition(A, &L, &U);
		}

		// matrix inversion
		{
			MatrixInverse(&iA, A);
		}

		// validation
		{
			MatrixMult(&ttA, A, iA);
		}

		printf("done!\n");
	}

	{
		// test newton's method
		ScalarF F = f_test_newton;

		Gradient g(2);
		Hessian H(2);

		g.Set(0, f_test_newton_g1);
		g.Set(1, f_test_newton_g2);

		H.Set(0, 0, f_test_newton_H00);
		H.Set(0, 1, f_test_newton_H01);
		H.Set(1, 0, f_test_newton_H10);
		H.Set(1, 1, f_test_newton_H11);

		ScalarVector initGuess(2);
		ScalarVector result(2);
		initGuess.Set(0, 0.75f);
		initGuess.Set(1, -1.25f);

		NewtonsMethod(F, &g, &H, initGuess, &result);
		printf("done!\n");
	}

	return 0;
}