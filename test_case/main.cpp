#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "../common/common_shared.h"
#include "../linear_algebra/scalar_matrix.h"
#include "../line_search/line_search_subp.h"
#include "../unconstrained_optimization/newtons_method.h"
#include "../constrained_optimization/lagrange_multiplier.h"

#include <Eigen/Dense>

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

void efunc1d12(const EVector& input, EVector* output)
{
	xassert(input.rows() == output->rows());
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	(*output)(0) = -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1);
	(*output)(1) = 200 * (x2 - x1 * x1);
}

//------------------------------------------------------------------------------------//
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

void efunc2d12(const EVector& input, EVector* output)
{
	xassert(input.rows() == output->rows());
	xassert(input.rows() == 2);
	float x1 = input(0);
	float x2 = input(1);

	(*output)(0) = -sinf(x1);
	(*output)(1) = -cosf(x2);
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

//------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------//
int main()
{

	//{
	//	LineSearchParams params;
	//	params.fMin = -1.f;
	//	params.rho = 0.01f;
	//	params.sigma = 0.1f;
	//	params.tau1 = 9.f;
	//	params.tau2 = 0.1f;
	//	params.tau3 = 0.5f;

	//	ScalarF F = func1;
	//	Gradient g(2);
	//	{
	//		g.Set(0, func1d1);
	//		g.Set(1, func1d2);
	//	}

	//	ScalarVector s(2); s.Set(0, 1.f); s.Set(1, 0.f);
	//	ScalarVector x0(2); x0.Set(0, 0.f); x0.Set(1, 0.f);

	//	float finalA = InexactLineSearch(F, g, s, x0, params);
	//	printf("done!\n");
	//}

	{
		LineSearchParams params;
		params.fMin = -1.f;
		params.rho = 0.01f;
		params.sigma = 0.1f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		const ScalarFunc F = efunc1;
		const GradientFunc g = efunc1d12;

		EVector s(2); s(0) = 1.f; s(1) = 0.f;
		EVector x0(2); x0(0) = 0.f; x0(1) = 0.f;

		float finalA = InexactLineSearch(F, g, s, x0, params);
		printf("done!\n");
	}

	//{
	//	LineSearchParams params;
	//	params.fMin = -2.f;
	//	params.rho = 0.01f;
	//	params.sigma = 0.001f;
	//	params.tau1 = 9.f;
	//	params.tau2 = 0.1f;
	//	params.tau3 = 0.5f;

	//	ScalarF F = func2;
	//	Gradient g(2);
	//	{
	//		g.Set(0, func2d1);
	//		g.Set(1, func2d2);
	//	}

	//	ScalarVector s(2); s.Set(0, 0.707f); s.Set(1, 0.707f);
	//	ScalarVector x0(2); x0.Set(0, 0.f); x0.Set(1, 0.f);

	//	float finalA = InexactLineSearch(F, g, s, x0, params);
	//	printf("done!\n");
	//}

	{
		LineSearchParams params;
		params.fMin = -2.f;
		params.rho = 0.01f;
		params.sigma = 0.001f;
		params.tau1 = 9.f;
		params.tau2 = 0.1f;
		params.tau3 = 0.5f;

		ScalarFunc F = efunc2;
		GradientFunc g = efunc2d12;

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

	//{
	//	// test newton's method
	//	ScalarF F = f_test_newton;

	//	Gradient g(2);
	//	Hessian H(2);

	//	g.Set(0, f_test_newton_g1);
	//	g.Set(1, f_test_newton_g2);

	//	H.Set(0, 0, f_test_newton_H00);
	//	H.Set(0, 1, f_test_newton_H01);
	//	H.Set(1, 0, f_test_newton_H10);
	//	H.Set(1, 1, f_test_newton_H11);

	//	ScalarVector initGuess(2);
	//	ScalarVector result(2);
	//	initGuess.Set(0, 0.75f);
	//	initGuess.Set(1, -1.25f);

	//	NewtonsMethodParams params;

	//	NewtonsMethod(F, &g, &H, params, initGuess, &result);
	//	printf("done!\n");
	//}

	//{
	//	// test newton's method
	//	ScalarFunc F = ef_test_newton;

	//	EGradient g(2);
	//	EHessian H(2);

	//	g.Set(0, ef_test_newton_g1);
	//	g.Set(1, ef_test_newton_g2);

	//	H.Set(0, 0, ef_test_newton_H00);
	//	H.Set(0, 1, ef_test_newton_H01);
	//	H.Set(1, 0, ef_test_newton_H10);
	//	H.Set(1, 1, ef_test_newton_H11);

	//	EVector initGuess(2);
	//	EVector result(2);
	//	initGuess(0) = 0.75f;
	//	initGuess(1) = -1.25f;

	//	NewtonsMethodParams params;

	//	NewtonsMethod(F, &g, &H, params, initGuess, &result);
	//	printf("done!\n");
	//}


	{
		auto func2 = [](const EVector& input)->float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return cosf(x1) - sinf(x2);
		};

		auto func2d12 = [](const EVector& input, EVector* output)
		{
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = -sinf(x1);
			(*output)(1) = -cosf(x2);
		};

		auto func2h12 = [](const EVector& input, EMatrix* output)
		{
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = -cosf(x1);
			(*output)(0, 1) = 0.f;
			(*output)(1, 0) = 0.f;
			(*output)(1, 1) = sinf(x2);
		};

		// test newton's method.
		ScalarFunc F = func2;
		GradientFunc g = func2d12;
		HessianFunc H = func2h12;

		EVector x0(2);
		EVector xstar(2);
		x0(0) = (-30.f * 3.1415f / 180.f);
		x0(1) = (-30.f * 3.1415f / 180.f);

		NewtonsMethodParams params;
		params.m_maxIter = 20;

		NewtonsMethod(F, g, H, params, x0, &xstar);
		printf("done!\n");
	}

	//{
	//	// test Symmetric Rank One method.
	//	ScalarF F = func3;
	//	Gradient g(2);
	//	g.Set(0, func3d1);
	//	g.Set(1, func3d2);

	//	ScalarVector x0(2);
	//	ScalarVector xstar(2);
	//	x0.Set(0, 0.1f);
	//	x0.Set(1, 1.f);

	//	NewtonsMethodParams params;
	//	params.m_min = -1.f;

	//	QuasiNewtonSR1(F, &g, params, x0, &xstar);
	//	printf("done!\n");
	//}

	//{
	//	ScalarF F = func4;
	//	Gradient g(6);
	//	g.Set(0, func4d1);
	//	g.Set(1, func4d2);
	//	g.Set(2, func4d3);
	//	g.Set(3, func4d4);
	//	g.Set(4, func4d5);
	//	g.Set(5, func4d6);

	//	ScalarVector x0(6);
	//	ScalarVector xstar0(6);
	//	ScalarVector xstar1(6);
	//	ScalarVector xstar2(6);
	//	x0.Set(0, 1.f);
	//	x0.Set(1, 1.f);
	//	x0.Set(2, 1.f);
	//	x0.Set(3, 1.f);
	//	x0.Set(4, 1.f);
	//	x0.Set(5, 1.f);

	//	NewtonsMethodParams params;
	//	params.m_min = -1.f;

	//	//QuasiNewtonSR1(F, &g, params, x0, &xstar0);
	//	QuasiNewtonDFP(F, &g, params, x0, &xstar1);
	//	QuasiNewtonBFGS(F, &g, params, x0, &xstar2);

	//	float test1 = F(xstar1);
	//	float test2 = F(xstar2);
	//	printf("done!\n");
	//}

	{
		auto func5 = [](const EVector& input)->float {
			xassert(input.rows() == 3);
			float a = input(0);
			float b = input(1);
			float c = input(2);

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
		};

		auto func5d123 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 3);
			float a = input(0);
			float b = input(1);
			float c = input(2);

			const int length = sizeof(sample5x) / sizeof(sample5x[0]);

			float sum0 = 0.f;
			float sum1 = 0.f;
			float sum2 = 0.f;

			for (int i = 0; i < length; i++)
			{
				float x = sample5x[i];
				float sample = sample5y[i];
				float t = (a * sinf(x) + b * cosf(x) + c - sample);
				sum0 += 2 * t * sinf(x);
				sum1 += 2 * t * cosf(x);
				sum2 += 2 * t;
			}

			(*output)(0) = sum0;
			(*output)(1) = sum1;
			(*output)(2) = sum2;
		};

		ScalarFunc F = func5;
		GradientFunc g = func5d123;

		EVector x0(3);
		EVector xstar0(3);
		EVector xstar1(3);
		EVector xstar2(3);
		x0(0) = 0.1f;
		x0(1) = 1.f;
		x0(2) = 0.f;

		NewtonsMethodParams params;
		params.m_min = -1.f;
		params.m_epsilon = 0.01f;

		//QuasiNewtonSR1(F, &g, params, x0, &xstar0);

		QuasiNewtonDFP(F, g, params, x0, &xstar1);
		QuasiNewtonBFGS(F, g, params, x0, &xstar2);

		const float ftest = F(xstar2);

		printf("done!\n");
	}

	{
		LagrangeMultMethodParams params;
		params.m_maxIter = 20;

		auto func6 = [](const EVector& input)-> float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1*x1*x1*x1 + 3 * x1*x1*x2 + 2 * x1*x2*x2 + x2*x2*x2*x2;
		};

		auto func6d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 4 * x1*x1*x1 + 6 * x1*x2 + 2 * x2*x2;
			(*output)(1) = 3 * x1*x1 + 4 * x1*x2 + 4 * x2*x2*x2;
		};

		auto func6h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 12 * x1*x1 + 6 * x2;
			(*output)(0, 1) = 6 * x1 + 4 * x2;
			(*output)(1, 0) = 6 * x1 + 4 * x2;
			(*output)(1, 1) = 4 * x1 + 12 * x2*x2;
		};

		auto cfunc6 = [](const EVector& input)->float {
			xassert(input.rows() >= 2);
			float x1 = input(0);
			float x2 = input(1);

			return (x1*x1 + x2*x2 - 1.f) * 1.f;
		};

		auto cfunc6d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 2 * x1 * 1.f;
			(*output)(1) = 2 * x2 * 1.f;
		};

		auto cfunc6h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 2 * 1.f;
			(*output)(0, 1) = 0 * 1.f;
			(*output)(1, 0) = 0 * 1.f;
			(*output)(1, 1) = 2 * 1.f;
		};

		ScalarFunc F = func6;
		GradientFunc gF = func6d12;
		HessianFunc hF = func6h12;

		ScalarFunc c = cfunc6;
		GradientFunc gC = cfunc6d12;
		HessianFunc hC = cfunc6h12;

		EVector x1(2); 
		x1(0) = 0.707f; x1(1) = sqrtf(1.f - x1(0) * x1(0));
		x1(0) = 0.f; x1(1) = 1.f;
		EVector xstar(2);
		
		LagrangeMultMethodResult res = LagrangeMultMethod(F, gF, hF, c, gC, hC, params, x1, &xstar);

		const float cstar = c(xstar);
		const float fstar = F(xstar);

		printf("done!\n");
	}

	{
		LagrangeMultMethodParams params;
		params.m_maxIter = 20;

		auto func7 = [](const EVector& input)->float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1*x1 + x2*x2;
		};

		auto func7d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 2 * x1;
			(*output)(1) = 2 * x2;
		};

		auto func7h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 2;
			(*output)(0, 1) = 0;
			(*output)(1, 0) = 0;
			(*output)(1, 1) = 2;
		};

		auto cfunc7 = [](const EVector& input)->float {
			xassert(input.rows() >= 2);
			float x1 = input(0);
			float x2 = input(1);

			return 2 * x1 + x2 - 2;
		};

		auto cfunc7d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 2;
			(*output)(1) = 1;
		};

		auto cfunc7h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 3);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 0;
			(*output)(0, 1) = 0;
			(*output)(1, 0) = 0;
			(*output)(1, 1) = 0;
		};

		ScalarFunc F = func7;
		GradientFunc gF = func7d12;
		HessianFunc hF = func7h12;

		ScalarFunc c = cfunc7;
		GradientFunc gC = cfunc7d12;
		HessianFunc hC = cfunc7h12;

		EVector x1(2); x1(0) = 1.f; x1(1) = 0.f;
		EVector xstar(2);

		LagrangeMultMethod(F, gF, hF, c, gC, hC, params, x1, &xstar);
		printf("done!\n");
	}

	{
		// minimize cost function x1 * x2^2
		// s.t: x1^2 + x2^2 - 2 = 0

		auto func8 = [](const EVector& input) -> float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1 * x2 * x2;
		};

		auto func8d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = x2*x2;
			(*output)(1) = 2 * x1*x2;
		};

		auto func8h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 0;
			(*output)(0, 1) = 2 * x2;
			(*output)(1, 0) = 2 * x2;
			(*output)(1, 1) = 2 * x1;
		};

		auto cfunc8 = [](const EVector& input) -> float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1*x1 + x2*x2 - 2.f;
		};

		auto cfunc8d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 2 * x1;
			(*output)(1) = 2 * x2;
		};

		auto cfunc8h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 2;
			(*output)(0, 1) = 0;
			(*output)(1, 0) = 0;
			(*output)(1, 1) = 2;
		};

		ScalarFunc F = func8;
		GradientFunc gF = func8d12;
		HessianFunc hF = func8h12;

		ScalarFunc c = cfunc8;
		GradientFunc gC = cfunc8d12;
		HessianFunc hC = cfunc8h12;

		EVector x1(2); x1(0) = -1.5f; x1(1) = -1.6f;
		EVector xstar(2);

		LagrangeMultMethodParams params;
		params.m_maxIter = 20;
		params.m_lamda1 = 0.f;

		SQP1(F, gF, hF, c, gC, hC, params, x1, &xstar);

		printf("done!\n");
	}

	{
		// minimize cost function x1 * x2^2
		// s.t: x1^2 + x2^2 - 2 <= 0

		auto func8 = [](const EVector& input) -> float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1 * x2 * x2;
		};

		auto func8d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = x2*x2;
			(*output)(1) = 2 * x1*x2;
		};

		auto func8h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 0;
			(*output)(0, 1) = 2 * x2;
			(*output)(1, 0) = 2 * x2;
			(*output)(1, 1) = 2 * x1;
		};

		auto cfunc8 = [](const EVector& input) -> float {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			return x1*x1 + x2*x2 - 2.f;
		};

		auto cfunc8d12 = [](const EVector& input, EVector* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0) = 2 * x1;
			(*output)(1) = 2 * x2;
		};

		auto cfunc8h12 = [](const EVector& input, EMatrix* output) {
			xassert(input.rows() == 2);
			float x1 = input(0);
			float x2 = input(1);

			(*output)(0, 0) = 2;
			(*output)(0, 1) = 0;
			(*output)(1, 0) = 0;
			(*output)(1, 1) = 2;
		};

		ScalarFunc F = func8;
		GradientFunc gF = func8d12;
		HessianFunc hF = func8h12;

		ScalarFunc c = cfunc8;
		GradientFunc gC = cfunc8d12;
		HessianFunc hC = cfunc8h12;

		EVector x1(2); x1(0) = -1.5f; x1(1) = -1.6f;
		EVector xstar(2);

		LagrangeMultMethodParams params;
		params.m_maxIter = 20;
		params.m_lamda1 = 0.f;

		SQP2(F, gF, hF, c, gC, hC, params, x1, &xstar);

		printf("SQP2 done!\n");
	}

	return 0;
}
