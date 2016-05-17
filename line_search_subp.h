
struct Interval
{
	float a;
	float b;

	Interval() {}
	Interval(float _a, float _b) : a(_a), b(_b) {}

	bool IsIn(float c) const;
};

struct FParams
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


BracketRes bracketing(float (*fa)(float), float (*fadev)(float), float a0, float a1, const FParams& params);
float sectioning(float (*fa)(float), float (*fadev)(float), const Interval& _prevI, const FParams& params);