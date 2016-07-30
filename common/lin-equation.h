// to solve a single linear equation.
struct LinEquation
{
	enum ReturnType
	{
		kExists,	// 
		kAlways,	// equation always true no matter which x, ex: 0 * x = 0.
		kNever,		// equation can never be satisfied no matter which x, ex: 0 * x = 1
	};

	enum Operator
	{
		kEq,	// equal to
		kLT,	// Less than
		kLTE,	// less than and equal
		kGT,	// greater than
		kGTE,	// greater than and equal
	};

	struct Result
	{
		ReturnType returnType;
		Operator op;
		float number;
	};

	static Operator OppositeOp(Operator op);
	// solve A * x "op" b
	static Result Solve(const float a, const float b, Operator op);
};

