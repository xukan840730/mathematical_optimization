
#include "../common/common_shared.h"
#include "math.h"

#include "line_search_subp.h"

// Min / Max
static inline float Min(float a, float b)
{
	return (a < b) ? a : b;
}

static inline float Max(float a, float b)
{
	return (a > b) ? a : b;
}

bool Interval::IsIn(float c) const
{
	float fMin = Min(a, b);
	float fMax = Max(a, b);

	return c >= fMin && c <= fMax;
}

// Hermite interpolating cubic in interval [0, 1]
// c(z) = c0 * z^0 + c1 * z^1 + c2 * z^2 + c3 * z^3, where
// c0 = f(0), 
// c1 = f'(0), 
// c2 = 3(f(1) - f(0)) - 2*f'(0) - f'(1)
// c3 = f'(0) + f'(1) - 2*(f(1) - f(0))

// alpha = a + z(b - a), df/dz = df/dalpha * dalpha/dz = (b - a) * df/dalpha

struct HermiteRes
{
	float c0;
	float c1;
	float c2;
	float c3;
};

// mapping from interval [a, b] to [0, 1]
HermiteRes hermite_interpolation01(const Interval& I, float fa, float fb, float devA, float devB)
{
	HermiteRes res;

	float t = (I.b - I.a);
	float fdev0 = t * devA;
	float fdev1 = t * devB;

	res.c0 = fa;
	res.c1 = fdev0;
	res.c2 = 3 * (fb - fa) - 2 * fdev0 - fdev1;
	res.c3 = fdev0 + fdev1 - 2 * (fb - fa);

	return res;
}

float choose(const Interval& currI, float currFa, float currFb, float currDevA, float currDevB,
	const Interval& nextI)
{
	HermiteRes hermiteRes = hermite_interpolation01(currI, currFa, currFb, currDevA, currDevB);

	// map alpha to z space.
	// the old interval was mapped to [0, 1], calculate the next interval's mapping.
	Interval mappingI;
	mappingI.a = (nextI.a - currI.a) / (currI.b - currI.a);
	mappingI.b = (nextI.b - currI.a) / (currI.b - currI.a);

	// calculate c(z) at each end of interval
	float cMappingA = hermiteRes.c0 + hermiteRes.c1 * mappingI.a + hermiteRes.c2 * mappingI.a * mappingI.a + hermiteRes.c3 * mappingI.a * mappingI.a * mappingI.a;
	float cMappingB = hermiteRes.c0 + hermiteRes.c1 * mappingI.b + hermiteRes.c2 * mappingI.b * mappingI.b + hermiteRes.c3 * mappingI.b * mappingI.b * mappingI.b;

	// d'z = c1 + 2*c2*z + 3*c3*z^2
	// d''z = 2*c2 + 6*c3*z, for a stationary point, need df/dz = 0, and d''z > 0.
	// so that (z + c2/(3*c3))^2 = (c2/c3)^2/9 - c1/(3*c3) ==> 

	float co = (hermiteRes.c2/hermiteRes.c3)*(hermiteRes.c2/hermiteRes.c3) / 9.f - hermiteRes.c1 / (3 * hermiteRes.c3);
	float selectedZ;

	if (co < 0.f)
	{
		selectedZ = cMappingA < cMappingB ? mappingI.a : mappingI.b;
	}
	else
	{
		float cZ0 = sqrtf(co) - hermiteRes.c2 / (3 * hermiteRes.c3);
		float cZ1 = -sqrtf(co) - hermiteRes.c2 / (3 * hermiteRes.c3);

		// second order derivative.
		float cZ0_dev2 = 2 * hermiteRes.c2 + 6 * hermiteRes.c3 * cZ0;
		float cZ1_dev2 = 2 * hermiteRes.c2 + 6 * hermiteRes.c3 * cZ1;

		if (cZ0_dev2 >= 0.f && mappingI.IsIn(cZ0))
		{
			selectedZ = cZ0;
		}
		else if (cZ1_dev2 >= 0.f && mappingI.IsIn(cZ1))
		{
			selectedZ = cZ1;
		}
		else
		{
			selectedZ = cMappingA < cMappingB ? mappingI.a : mappingI.b;
		}
	}

	// unmap z to alpha space.
	float selectedAlpha = currI.a + (currI.b - currI.a) * selectedZ;
	return selectedAlpha;
}


BracketRes bracketing(float (*fa)(float), float (*fadev)(float), float a0, float a1, const LineSearchParams& params)
{
	// a0 == 0.
	xassert(a0 == 0);

	BracketRes res;

	float Ai = a1;

	float FA0 = fa(a0);
	float FAD0 = fadev(a0);

	float AiMinus1 = a0;
	float FAiMinus1 = FA0;
	float FADevIMinus1 = FAD0;

	xassert(FAD0 <= 0.f);

	float mu = (params.fMin - FA0) / (params.rho * FAD0);
	xassert(mu > a0);

	while (true)
	{
		// evaluate f(Ai)
		float FAi = fa(Ai);

		if (FAi <= params.fMin)
		{
			res.t = true;
			res.alpha = Ai;
			return res;
		}

		{
			bool tB0 = FAi > FA0 + Ai * params.rho * FAD0;
			bool tB1 = FAi >= FAiMinus1;
			if (tB0 || tB1)
			{
				res.tB = true;
				res.interval = Interval(AiMinus1, Ai);
				return res;
			}
		}

		// evalulate f'(Ai)
		float FADevI = fadev(Ai);

		{
			if (fabs(FADevI) <= -params.sigma * FAD0)
			{
				res.t = true;
				res.alpha = Ai;
				return res;
			}

			if (FADevI >= 0)
			{
				res.tB = true;
				res.interval = Interval(Ai, AiMinus1);
				return res;
			}
		}

		{
			float AiPlus1 = Ai;

			if (mu <= 2 * Ai - AiMinus1)
			{
				AiPlus1 = mu;
			}
			else
			{
				float M0 = 2 * Ai - AiMinus1;
				float M1 = Min(mu, Ai + params.tau1 * (Ai - AiMinus1));

				//float FM0 = fa(M0);
				//float FM1 = fa(M1);
				//float FMdev0 = fadev(M0);
				//float FMdev1 = fadev(M1);
				//AiPlus1 = choose(M0, M1, FM0, FM1, FMdev0, FMdev1);

				Interval currI;
				currI.a = AiMinus1;
				currI.b = Ai;

				Interval nextI;
				nextI.a = M0;
				nextI.b = M1;

				AiPlus1 = choose(currI, FAiMinus1, FAi, FADevIMinus1, FADevI, nextI);
			}

			// save last Ai, FAi and FDevAi
			AiMinus1 = Ai;
			FAiMinus1 = FAi;
			FADevIMinus1 = FADevI;

			// update Ai => A(i+1)
			Ai = AiPlus1;
		}
	}

	return res;
}

float sectioning(float (*fa)(float), float (*fadev)(float), const Interval& _prevI, const LineSearchParams& params)
{
	xassert(_prevI.a == _prevI.a);
	xassert(_prevI.b == _prevI.b);

	Interval currJ = _prevI;

	const float FA0 = fa(0.f);
	const float FAD0 = fadev(0.f);

	while (true)
	{
		const float dJ = currJ.b - currJ.a;
		xassert(dJ != 0.f);

		// choose alphaj belongs to [aj + tau2(bj - aj), bj - tau3(bj - aj)];
		Interval range;
		{
			range.a = currJ.a + params.tau2 * dJ;
			range.b = currJ.b - params.tau3 * dJ;
		}

		float alphaJ = choose(currJ, fa(currJ.a), fa(currJ.b), fadev(currJ.a), fadev(currJ.b), range);
		xassert(alphaJ == alphaJ);

		// evaluate f(alphaJ)
		float FAj = fa(alphaJ);

		bool t0 = FAj > FA0 + params.rho * alphaJ * FAD0;
		bool t1 = FAj >= fa(currJ.a);
		if (t0 || t1)
		{
			currJ.a = currJ.a;
			currJ.b = alphaJ;
		}
		else
		{
			float FADj = fadev(alphaJ);
			if (fabs(FADj) <= -params.sigma * FAD0)
			{
				// finally we get a good point.
				return alphaJ;
			}
			else
			{
				if (dJ * FADj >= 0)
				{
					float oa = currJ.a;
					currJ.a = alphaJ;
					currJ.b = oa;
				}
				else
				{
					currJ.a = alphaJ;
					currJ.b = currJ.b;
				}
			}
		}
	}

	return 0.f;
}
