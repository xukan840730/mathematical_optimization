#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "qpsub.h"
#include "eqnsolv.h"

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

void qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector& lb, const EVector& ub, const EVector* x0, int numEqCstr, QpsubCaller caller, float eps) 
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
			A = MatRowAppend(A, MatFromRowIdx(lbmat, arglb));
			b = VecAppend(b, VecFromIdx(-lb, arglb));
		}
	}
	// upper bound
	const EVector argub = findUb(ub);
	{
		if (argub.rows() > 0)
		{
			EMatrix ubmat(numVars, numVars); ubmat.setIdentity(); // ubmat = I
			A = MatRowAppend(A, MatFromRowIdx(ubmat, argub));
			b = VecAppend(b, VecFromIdx(ub, argub));
		}
	}
	numCstr = numCstr + arglb.rows() + argub.rows();
	ASSERT(A.rows() == numCstr);

	int maxIter = 200; // TODO:
	// used for determining threshold for whether a direction will violate
	EVector normA(numCstr); normA.setOnes();
	if (normalize)
	{
		for (int ii = 0; ii < numCstr; ii++)
		{
			float rowNorm = A.row(ii).norm();
			if (rowNorm > 0.f)
			{
				A.row(ii) /= rowNorm;
				b(ii) /= rowNorm;
				normA(ii) = rowNorm;
			}
		}
	}
	// some error number
	float errNorm = 0.01f * sqrt(eps);
	float tolDep = 100 * numVars * eps;
	EVector lambda(numCstr); lambda.setZero();
	EVector eqix = colon(0, numEqCstr - 1);
	EVector indepIdx = colon(0, numCstr - 1); // independent constraint indices

	if (numEqCstr > 0)
	{
		eqnres res = eqnsolv(A, b, eqix, numVars, eps); 
		if (res.rmvIdx.rows() > 0)
		{
			indepIdx = VecRmvIdx(indepIdx, res.rmvIdx);
			normA = VecFromIdx(normA, indepIdx);
		}
	}
	else
	{}
}
