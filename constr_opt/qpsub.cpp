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

int qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector& lb, const EVector& ub, const EVector* _x0, int numEqCstr, QpsubCaller caller, float eps) 
{
	ASSERT(_A.rows() == _b.rows());
	if (_x0 != nullptr) ASSERT(_x0->rows() == H.cols());
	int exitFlag = 1;
	int iteration = 0;

	// copy original matrix and vector.
	EVector f = _f;
	EMatrix A = _A;
	EVector b = _b;

	int numCstr = A.rows();
	
	bool lls = false;
	bool normalize = false;
	int numVars = _x0 != nullptr ? _x0->rows() : H.cols();
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

	int maxIter = 200; // TODO: pass in as a parameter
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
	const float errNorm = 0.01f * sqrt(eps);
	const float tolDep = 100 * numVars * eps;
	const float tolCons = 1e-10;
	EVector lambda(numCstr); lambda.setZero();
	EVector eqix = colon(0, numEqCstr - 1);
	EVector indepIdx = colon(0, numCstr - 1); // independent constraint indices

	EVector X0;
	if (_x0 != nullptr) { X0 = *_x0; }
	else { X0 = EVector(numVars); X0.setZero(); }

	int actCnt = 0;
	if (numEqCstr > 0)
	{
		eqnres res = eqnsolv(A, b, eqix, numVars, eps);
		if (res.rmvIdx.rows() > 0)
		{
			indepIdx = VecRmvIdx(indepIdx, res.rmvIdx);
			normA = VecFromIdx(normA, indepIdx);
		}

		if (res.exitFlag == -1)
		{
			// equalities are inconsistent, so x and lamabda have no valid values
			// return original x and zero for lambda
			return res.exitFlag;
		}

		// update flags.
		numEqCstr = eqix.rows();
		numCstr = A.rows();

		const EMatrix Aeq = MatFromRowIdx(A, eqix);
		const EVector beq = VecFromIdx(b, eqix);
		// is this necessary?
		{
			// find a feasible point for equality constraints
			EVector diff = Aeq * X0 - beq;
			if (VecAbs(diff).maxCoeff() > tolCons)
				X0 = EigenColPivQrSolve(Aeq, beq);
		}

		// is this necessary?
		if (numEqCstr > numVars)
		{
			EVector diff = Aeq * X0 - beq;
			float err = VecAbs(diff).maxCoeff();
			if (err > 1e-8) // equalities not met
			{
				// exiting: the equality constraints are overly stringent
				// there's no feasible solution
				exitFlag = -1;
				return exitFlag;
			}
			else
			{
				// check inequalities
				EVector diff2 = A * X0 - b;
				if (diff2.maxCoeff() > 1e-8)
				{
					// exiting: the constraints or bounds are overly stringent
					// there's no feasible solution
					// equality constraints have been met
					exitFlag = -1;
				}
			}
		}

		EMatrix Q, R;
		EigenQrDecomp(Aeq, &Q, &R, nullptr);

		EVector actLambda;
		if (isqp)
			actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H*X0 + f));
		else if (lls)
			actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H.transpose() * (H*X0 - f)));
		else
			actLambda = EigenColPivQrSolve(-R, Q.transpose() * f);
		
		lambda(indepInd(eqix)) = normf * (actLambda ./ normA(eqix));
		return exitFlag;
	}
	else
	{}

	return 0;
}
