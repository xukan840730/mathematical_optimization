#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "qpsub.h"

static bool validLb(float x, void* params)
{
	float t = *(static_cast<float*>(params));
	return t > -NDI_FLT_MAX;
}

static EVector findLb(const EVector& lb)
{
	float t = -NDI_FLT_MAX;
	return findRows(lb, validLb, &t);
}

void qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector& lb, const EVector& ub, const EVector* x0, int numEqCstr, int numCstr, QpsubCaller caller) 
{
	ASSERT(_A.rows() == _b.rows());
	if (x0 != nullptr) ASSERT(x0->rows() == H.cols());
	int exitFlag = 1;
	int iteration = 0;

	// copy original matrix and vector.
	EVector f = _f;
	EMatrix A = _A;
	EVector b = _b;
	
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
	const EVector arglb = findLb(lb);
	if (anyNnz(arglb, NDI_FLT_EPSILON))
	{
		EMatrix lbmat(numVars, numVars); lbmat.setIdentity();
		lbmat *= -1;
		A = MatRowAppend(A, MatrixFromRowIdx(lbmat, arglb));
		b = VecAppend(b, VectorFromIdx(-lb, arglb));
	}
}
