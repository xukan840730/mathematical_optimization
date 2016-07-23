#ifndef _QUADRATIC_PROGRAMMING_H_
#define _QUADRATIC_PROGRAMMING_H_

struct FeasibilityRes
{
public:
	enum Type {
		kFeasible,		// the minimal of EQP is feasible in original QP problem.
		kInfeasible,	// the minimal of EQP violates constraints of original QP problem
		kImpossible,	// it's impossible to find t between [0-1] to move along x0 to xeqp
	};

	Type type;
	float t;
	int vconstrIdx;		// the first violated constrained index.
};

// to find a feasible t from x0 to xp and satisfy linear constraints Ax <= b.
FeasibilityRes LConstrFeasibility(const EVector& xstart, const EVector& xend, const EMatrix& A, const EVector& b);

struct EQuadProgRes
{
	bool success;
	EVector xstar;
};
// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. A.x = b
EQuadProgRes EQuadProg(const EMatrix& H, const EVector& q, const EMatrix& Aeq, const EVector& beq);

// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. A.x <= b
void QuadProg(const EMatrix& H, const EVector& q, const EMatrix& A, const EVector& b);
// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. (1) A.x <= b, Aeq.x = beq.
void QuadProg(const EMatrix& H, const EVector& q, const EMatrix& A, const EVector& b, const EMatrix& Aeq, const EVector& beq);

#endif
