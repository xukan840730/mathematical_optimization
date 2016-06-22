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

//------------------------------------------------------------------------------------//
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

//------------------------------------------------------------------------------------//
float func1(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
}

float func1d1(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1);
}

float func1d2(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 200 * (x2 - x1 * x1);
}

//------------------------------------------------------------------------------------//
float func2(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);
	
	return cosf(x1) - sinf(x2);
}

float func2d1(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return -sinf(x1);
}

float func2d2(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return -cosf(x2);
}

float func2h00(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return -cosf(x1);
}

float func2h01(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 0.f;
}

float func2h10(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 0.f;
}

float func2h11(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return sinf(x2);
}

//------------------------------------------------------------------------------------//
float func3(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 10 * x1 * x1 + x2 * x2;
}

float func3d1(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 20 * x1;
}

float func3d2(const ScalarVector& input)
{
	xassert(input.GetLength() == 2);
	float x1 = input.Get(0);
	float x2 = input.Get(1);

	return 2 * x2;
}

//------------------------------------------------------------------------------------//
<<<<<<< HEAD
float sample4x[4] = {-1.f, 1.f, 2.f, 4.f};
float sample4y[4] = {1.f, -1.f, 1.f, 6.f};

// f() = a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
float func4(const ScalarVector& input)
{
	// from f() to get cost func.
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += t * t;
	}

	return sum;
}

float func4d1(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t * (x * x * x * x * x);
	}

	return sum;
}

float func4d2(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t * (x * x * x * x);
	}

	return sum;
}

float func4d3(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t * (x * x * x);
	}

	return sum;
}

float func4d4(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t * (x * x);
	}

	return sum;
}

float func4d5(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t * (x);
	}

	return sum;
}

float func4d6(const ScalarVector& input)
{
	xassert(input.GetLength() == 6);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);
	float d = input.Get(3);
	float e = input.Get(4);
	float f = input.Get(5);

	float sum = 0.f;

	const int length = sizeof(sample4x) / sizeof(sample4x[0]);
	for (int i = 0; i < length; i++)
	{
		float x = sample4x[i];
		float sample = sample4y[i];
		float t = (a * x * x * x * x * x + b * x * x * x * x + c * x * x * x + d * x * x + e * x + f - sample);
		sum += 2 * t;
	}

=======
float func4x[4] = {-1.f, 1.f, 2.f, 4.f};
float func4y[4] = {1.f, -1.f, 1.f, 6.f};

float func4(const ScalarVector& input)
{
	// f(a, b, c) = a*x^2  + b*x + c
	// cost func.
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	float sum = 0.f;
	for (int i = 0; i < 4; i++)
	{
		float x = func4x[i];
		float t = (a * x * x + b * x + c - func4y[i]);
		sum += t * t;
	}
	return sum;
}

float func4da(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	float sum = 0.f;
	for (int i = 0; i < 4; i++)
	{
		float x = func4x[i];
		float t = (a * x * x + b * x + c - func4y[i]);
		sum += 2 * t * x * x;
	}
	return sum;
}

float func4db(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	float sum = 0.f;
	for (int i = 0; i < 4; i++)
	{
		float x = func4x[i];
		float t = (a * x * x + b * x + c - func4y[i]);
		sum += 2 * t * x;
	}
	return sum;
}

float func4dc(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	float sum = 0.f;
	for (int i = 0; i < 4; i++)
	{
		float x = func4x[i];
		float t = (a * x * x + b * x + c - func4y[i]);
		sum += 2 * t;
	}
>>>>>>> c8a15478aa90351f14259197d47fe8953b2d92c7
	return sum;
}


//------------------------------------------------------------------------------------//
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

		ScalarF F = func1;
		Gradient g(2);
		{
			g.Set(0, func1d1);
			g.Set(1, func1d2);
		}

		ScalarVector s(2); s.Set(0, 1.f); s.Set(1, 0.f);
		ScalarVector x0(2); x0.Set(0, 0.f); x0.Set(1, 0.f);

		float finalA = InexactLineSearch(F, g, s, x0, params);
		printf("done!\n");
	}

	{
		LineSearchParams params;
		params.fMin = -2.f;
		params.rho = 0.01f;
		params.sigma = 0.001f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		ScalarF F = func2;
		Gradient g(2);
		{
			g.Set(0, func2d1);
			g.Set(1, func2d2);
		}

		ScalarVector s(2); s.Set(0, 0.707f); s.Set(1, 0.707f);
		ScalarVector x0(2); x0.Set(0, 0.f); x0.Set(1, 0.f);

		float finalA = InexactLineSearch(F, g, s, x0, params);
		printf("done!\n");
	}

	{
		ScalarMatrix A(3, 3);
		ScalarMatrix iA(3, 3);

		ScalarMatrix ttA(3, 3);

		ScalarMatrix C(3, 3);
		ScalarMatrix LC(3, 3);
		ScalarMatrix tLC(3, 3);
		ScalarMatrix ttC(3, 3);

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

		// Cholesky decomposition
		{
			C.Set(0, 0, 6.f);
			C.Set(0, 1, 15.f);
			C.Set(0, 2, 55.f);

			C.Set(1, 0, 15.f);
			C.Set(1, 1, 55.f);
			C.Set(1, 2, 225.f);

			C.Set(2, 0, 55.f);
			C.Set(2, 1, 225.f);
			C.Set(2, 2, 979.f);

			CholeskyDecomposition(C, &LC);
			LC.Transpose(&tLC);

			MatrixMult(&ttC, LC, tLC);
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

		NewtonsMethodParams params;

		NewtonsMethod(F, &g, &H, params, initGuess, &result);
		printf("done!\n");
	}

	{
		// test newton's method.
		ScalarF F = func2;

		Gradient g(2);
		Hessian H(2);

		g.Set(0, func2d1);
		g.Set(1, func2d2);

		H.Set(0, 0, func2h00);
		H.Set(0, 1, func2h01);
		H.Set(1, 0, func2h10);
		H.Set(1, 1, func2h11);

		ScalarVector x0(2);
		ScalarVector xstar(2);
		x0.Set(0, -30.f * 3.1415f / 180.f);
		x0.Set(1, -30.f * 3.1415f / 180.f);

		NewtonsMethodParams params;
		params.m_maxIter = 20;

		NewtonsMethod(F, &g, &H, params, x0, &xstar);
		printf("done!\n");
	}

	{
		// test Symmetric Rank One method.
		ScalarF F = func3;
		Gradient g(2);
		g.Set(0, func3d1);
		g.Set(1, func3d2);

		ScalarVector x0(2);
		ScalarVector xstar(2);
		x0.Set(0, 0.1f);
		x0.Set(1, 1.f);

		NewtonsMethodParams params;
		params.m_min = -1.f;

		QuasiNewtonSR1(F, &g, params, x0, &xstar);
		printf("done!\n");
	}

	{
<<<<<<< HEAD
		ScalarF F = func4;
		Gradient g(6);
		g.Set(0, func4d1);
		g.Set(1, func4d2);
		g.Set(2, func4d3);
		g.Set(3, func4d4);
		g.Set(4, func4d5);
		g.Set(5, func4d6);

		ScalarVector x0(6);
		ScalarVector xstar0(6);
		ScalarVector xstar1(6);
		x0.Set(0, 1.f);
		x0.Set(1, 1.f);
		x0.Set(2, 1.f);
		x0.Set(3, 1.f);
		x0.Set(4, 1.f);
		x0.Set(5, 1.f);

		NewtonsMethodParams params;
		params.m_min = -1.f;

		QuasiNewtonSR1(F, &g, params, x0, &xstar0);
		QuasiNewtonDFP(F, &g, params, x0, &xstar1);
=======
		// use QuasiNewton method to test curve fitting
		ScalarF F = func4;
		Gradient g(3);
		g.Set(0, func4da);
		g.Set(1, func4db);
		g.Set(2, func4dc);

		ScalarVector x0(3);
		ScalarVector xstar(3);

		x0.Set(0, 1);
		x0.Set(1, 1);
		x0.Set(2, 0);

		NewtonsMethodParams params;
		params.m_min = NDI_FLT_EPSILON;

		QuasiNewtonSR1(F, &g, params, x0, &xstar);

		float fmin = F(xstar);

		ScalarVector xstar0 = xstar;
		ScalarVector xstar1 = xstar;
		xstar0.Multiply(0.999f);
		xstar1.Multiply(1.001f);

		float fmin0 = F(xstar0);
		float fmin1 = F(xstar1);
>>>>>>> c8a15478aa90351f14259197d47fe8953b2d92c7
		printf("done!\n");
	}

	return 0;
}