#ifndef _QUADRATIC_PROGRAMMING_H_
#define _QUADRATIC_PROGRAMMING_H_

struct EQuadProgRes
{
	enum {
		kFound,
		kUnbounded,
	};

	int type;
	EVector xstar;	// if minima found
	EVector stepp;	// if unbounded, stepp is the most negative curvature direction
};
// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. A.x = b
EQuadProgRes EQuadProg(const EMatrix& H, const EVector& q, const EMatrix& Aeq, const EVector& beq);

// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. Ain.x <= bin
int QuadProg(const EMatrix& H, const EVector& q, const EMatrix& Ain, const EVector& bin);
// minimize a quadratic problem: min 1/2 x^t H x + q^t x, s.t. (1) Ain.x <= bin, Aeq.x = beq.
//void QuadProg(const EMatrix& H, const EVector& q, const EMatrix& Ain, const EVector& bin, const EMatrix& Aeq, const EVector& beq);

#endif
