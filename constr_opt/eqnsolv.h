#ifndef _EQNSOLV_H_
#define _EQNSOLV_H_

// check equality constraints dependency and consistency.
// return -1 if eqcstr are inconsistent, or, remove duplicated constraints
struct eqnres
{
	int exitFlag;
	EVector rmvIdx;
};
eqnres eqnsolv(EMatrix& A, EVector& b, EVector& eqix, int numVars, float eps);

#endif
