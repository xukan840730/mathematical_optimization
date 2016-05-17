#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "line_search_subp.h"

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
	return cos(a) - sin(a);
}

float h_alpha_dev(float a)
{
	return -sin(a) - cos(a);
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


	return 0;
}