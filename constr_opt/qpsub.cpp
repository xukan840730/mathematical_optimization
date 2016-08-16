#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "qpsub.h"

static bool validLb(float x, void* params)
{
	float t = *(static_cast<float*>(params));
	return x > -NDI_FLT_MAX;
}

static EVector findLb(const EVector& lb)
{
	float t = -NDI_FLT_MAX;
	return findRows(lb, validLb, &t);
}

static bool validUb(float x, void* params)
{
	float t = *(static_cast<float*>(params));
	return x < NDI_FLT_MAX;
}

static EVector findUb(const EVector& ub)
{
	float t = NDI_FLT_MAX;
	return findRows(ub, validUb, &t);
}

void qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector& lb, const EVector& ub, const EVector* x0, int numEqCstr, QpsubCaller caller) 
{
	ASSERT(_A.rows() == _b.rows());
	if (x0 != nullptr) ASSERT(x0->rows() == H.cols());
	int exitFlag = 1;
	int iteration = 0;

	// copy original matrix and vector.
	EVector f = _f;
	EMatrix A = _A;
	EVector b = _b;

	int numCstr = A.rows();
	
	bool lls = false;
	bool normalize = false;
	int numVars = x0 != nullptr ? x0->rows() : H.cols();
	int simplexIter = 0;

	bool isqp = H.squaredNorm() > 0.f;
	if (caller == kLsqlin)
	{
		lls = true;
		isqp = false;
		normalize = true;
	}
	else if (caller == kQpsub)
	{
		normalize = false;
	}

	float normf = 1.f;
	if (normalize && !isqp && !lls)
	{
		normf = f.norm();
		if (normf > 0.f)
			f = f / normf;
	}

	// Handle bounds as linear constraints
	ASSERT(lb.rows() == numVars);
	ASSERT(ub.rows() == numVars);
	// lower bound
 	const EVector arglb = findLb(lb);
	{
		if (arglb.rows() > 0)
		{
			EMatrix lbmat(numVars, numVars); lbmat.setIdentity(); lbmat *= -1; // lbmat = -I
			A = MatRowAppend(A, MatrixFromRowIdx(lbmat, arglb));
			b = VecAppend(b, VectorFromIdx(-lb, arglb));
		}
	}
	// upper bound
	const EVector argub = findUb(ub);
	{
		if (argub.rows() > 0)
		{
			EMatrix ubmat(numVars, numVars); ubmat.setIdentity(); // ubmat = I
			A = MatRowAppend(A, MatrixFromRowIdx(ubmat, argub));
			b = VecAppend(b, VectorFromIdx(ub, argub));
		}
	}
	numCstr = numCstr + arglb.rows() + argub.rows();
}
