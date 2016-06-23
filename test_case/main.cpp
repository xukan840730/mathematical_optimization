#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "../common/common_shared.h"
#include "../linear_algebra/scalar_matrix.h"
#include "../line_search/line_search_subp.h"
#include "../newtons_method/newtons_method.h"

#include <Eigen/Dense>

using Eigen::MatrixXd;

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

float ef_test_newton(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return x1 * x1 * x1 * x1 + x1 * x2 + (1 + x2) * (1 + x2);
}

float ef_test_newton_g1(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return 4 * x1 * x1 * x1 + x2;
}

float ef_test_newton_g2(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return x1 + 2 * (1 + x2);
}

float ef_test_newton_H00(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return 12 * x1 * x1;
}

float ef_test_newton_H01(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return 1;
}

float ef_test_newton_H10(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
	return 1;
}

float ef_test_newton_H11(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);
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

float efunc1(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
}

float efunc1d1(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1);
}

float efunc1d2(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

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

float efunc2(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return cosf(x1) - sinf(x2);
}

float efunc2d1(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return -sinf(x1);
}

float efunc2d2(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return -cosf(x2);
}

float efunc2h00(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return -cosf(x1);
}

float efunc2h01(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return 0.f;
}

float efunc2h10(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	return 0.f;
}

float efunc2h11(const EVector& input)
{
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

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

	return sum;
}

//------------------------------------------------------------------------------------//
float sample5x[6] = {-1.f, 0.f, 1.f, 2.f, 3.f, 4.f};
float sample5y[6] = {1.f, 0.f, -1.f, 1.f, 2.5f, 6.f};

float func5(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	const int length = sizeof(sample5x) / sizeof(sample5x[0]);

	float sum = 0.f;
	for (int i = 0; i < length; i++)
	{
		float x = sample5x[i];
		float sample = sample5y[i];
		float t = (a * sinf(x) + b * cosf(x) + c - sample);
		sum += t * t;
	}
	return sum;
}

float func5d0(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	const int length = sizeof(sample5x) / sizeof(sample5x[0]);

	float sum = 0.f;
	for (int i = 0; i < length; i++)
	{
		float x = sample5x[i];
		float sample = sample5y[i];
		float t = (a * sinf(x) + b * cosf(x) + c - sample);
		sum += 2 * t * sinf(x);
	}
	return sum;
}

float func5d1(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	const int length = sizeof(sample5x) / sizeof(sample5x[0]);

	float sum = 0.f;
	for (int i = 0; i < length; i++)
	{
		float x = sample5x[i];
		float sample = sample5y[i];
		float t = (a * sinf(x) + b * cosf(x) + c - sample);
		sum += 2 * t * cosf(x);
	}
	return sum;
}

float func5d2(const ScalarVector& input)
{
	xassert(input.GetLength() == 3);
	float a = input.Get(0);
	float b = input.Get(1);
	float c = input.Get(2);

	const int length = sizeof(sample5x) / sizeof(sample5x[0]);

	float sum = 0.f;
	for (int i = 0; i < length; i++)
	{
		float x = sample5x[i];
		float sample = sample5y[i];
		float t = (a * sinf(x) + b * cosf(x) + c - sample);
		sum += 2 * t;
	}
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
		params.fMin = -1.f;
		params.rho = 0.01f;
		params.sigma = 0.1f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		ScalarFunc F = efunc1;
		EGradient g(2);
		{
			g.Set(0, efunc1d1);
			g.Set(1, efunc1d2);
		}

		EVector s(2); s(0) = 1.f; s(1) = 0.f;
		EVector x0(2); x0(0) = 0.f; x0(1) = 0.f;

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
		LineSearchParams params;
		params.fMin = -2.f;
		params.rho = 0.01f;
		params.sigma = 0.001f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		ScalarFunc F = efunc2;
		EGradient g(2);
		{
			g.Set(0, efunc2d1);
			g.Set(1, efunc2d2);
		}

		EVector s(2); s(0) = 0.707f; s(1) = 0.707f;
		EVector x0(2); x0(0) = 0.f; x0(1) = 0.f;

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
		// test newton's method
		ScalarFunc F = ef_test_newton;

		EGradient g(2);
		EHessian H(2);

		g.Set(0, ef_test_newton_g1);
		g.Set(1, ef_test_newton_g2);

		H.Set(0, 0, ef_test_newton_H00);
		H.Set(0, 1, ef_test_newton_H01);
		H.Set(1, 0, ef_test_newton_H10);
		H.Set(1, 1, ef_test_newton_H11);

		EVector initGuess(2);
		EVector result(2);
		initGuess(0) = 0.75f;
		initGuess(1) = -1.25f;

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
		ScalarVector xstar2(6);
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
		QuasiNewtonBFGS(F, &g, params, x0, &xstar2);
		printf("done!\n");
	}

	{
		ScalarF F = func5;
		Gradient g(3);
		g.Set(0, func5d0);
		g.Set(1, func5d1);
		g.Set(2, func5d2);

		ScalarVector x0(3);
		ScalarVector xstar0(3);
		ScalarVector xstar1(3);
		ScalarVector xstar2(3);
		x0.Set(0, 0.1f);
		x0.Set(1, 1.f);
		x0.Set(2, 0.f);

		NewtonsMethodParams params;
		params.m_min = -1.f;

		//QuasiNewtonSR1(F, &g, params, x0, &xstar0);
		QuasiNewtonDFP(F, &g, params, x0, &xstar1);
		QuasiNewtonBFGS(F, &g, params, x0, &xstar2);

		const float ftest = F(xstar2);

		printf("done!\n");
	}

	return 0;
}