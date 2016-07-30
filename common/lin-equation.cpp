#include "lin-equation.h"
#include <math.h>

LinEquation::Operator LinEquation::OppositeOp(Operator op)
{	
	switch (op)
	{
	case kEq: return kEq; 
	case kLT: return kGT;
	case kLTE: return kGTE;
	case kGT: return kLT;
	case kGTE: return kLTE;
	}
	return kEq;
}

LinEquation::Result LinEquation::Solve(const float a, const float b, Operator op)
{
	Result res;
	if (a == 0.f)
	{
		// first, if a is close to, or equal to 0.
		switch (op)
		{
		case kEq: res.returnType = b == 0.f ? kAlways : kNever; break;
		case kLT: res.returnType = b > 0.f ? kAlways : kNever; break; 
		case kLTE: res.returnType = b >= 0.f ? kAlways : kNever; break;
		case kGT: res.returnType = b < 0.f ? kAlways : kNever; break;
		case kGTE: res.returnType = b <= 0.f ? kAlways : kNever; break;
		}
	}
	else if (a < 0.f)
	{
		res.op = OppositeOp(op);
		res.number = b / a;
		res.returnType = kExists;
	}
	else 
	{
		res.op = op;
		res.number = b / a;
		res.returnType = kExists;
	}

	return res;
}
