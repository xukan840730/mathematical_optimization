
#include "../linear_algebra/gradient.h"

struct Interval
{
	float a;
	float b;

	Interval() {}
	Interval(float _a, float _b) : a(_a), b(_b) {}

	bool IsIn(float c) const;
};

struct LineSearchParams
{
	float fMin;
	float rho;
	float sigma;
	float tau1;
	float tau2;
	float tau3;
};

struct BracketRes
{
	BracketRes()
		: t(false)
		, tB(false)
	{}

	bool t;		// algorithm finishes with a selected alpha.
	bool tB;	// algorithm finishes with an interval.
	float alpha;

	Interval interval;
};

BracketRes Bracketing(const ScalarF F, const Gradient& g, const ScalarVector& s, 
	const ScalarVector& x0, float a1, const LineSearchParams& params);

BracketRes Bracketing(const ScalarFunc F, const GradientFunc g, const EVector& s, 
	const EVector& x0, float a1, const LineSearchParams& params);

float Sectioning(const ScalarF F, const Gradient& g, const ScalarVector& s,
	const ScalarVector& x0, const Interval& _prevI, const LineSearchParams& params);

float Sectioning(const ScalarFunc F, const GradientFunc g, const EVector& s,
	const EVector& x0, const Interval& _prevI, const LineSearchParams& params);

// inexact line search method. returns alpha along search direction s.
float InexactLineSearch(const ScalarF F, const Gradient& g, const ScalarVector& s, const ScalarVector& x0, const LineSearchParams& params);
float InexactLineSearch(const ScalarFunc F, const GradientFunc g, const EVector& s, const EVector& x0, const LineSearchParams& params);