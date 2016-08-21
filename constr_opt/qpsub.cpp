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

static void CalcInfeasibleLambda(bool isqp, bool lls, const EMatrix& Q, const EMatrix& R, const EMatrix& H, const EVector& f, const EVector& X, const EVector& eqix, const EVector& indepIdx, float normf, const EVector& normA, EVector& lambda)
{
	EVector actLambda;
	if (isqp)
		actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H*X + f));
	else if (lls)
		actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H.transpose() * (H*X - f)));
	else
		actLambda = EigenColPivQrSolve(-R, Q.transpose() * f);
	const EVector normedActLambda = normf * VecDivVec(actLambda, VecFromIdx(normA, eqix));
	VecChangeRows(lambda, normedActLambda, VecFromIdx(indepIdx, eqix));
}

int qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector* _lb, const EVector* _ub, const EVector* _x0, int numEqCstr, QpsubCaller caller, float eps) 
{
	ASSERT(_A.rows() == _b.rows());
	if (_x0 != nullptr) ASSERT(_x0->rows() == H.cols() || H.cols() == 0);
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
 	EVector arglb;
	if (_lb != nullptr)
	{
		// lower bound
		const EVector& lb = *_lb;
		ASSERT(lb.rows() == numVars);
		arglb = findLb(lb);
		if (arglb.rows() > 0)
		{
			EMatrix lbmat(numVars, numVars); lbmat.setIdentity(); lbmat *= -1; // lbmat = -I
			A = MatRowAppend(A, MatFromRowIdx(lbmat, arglb));
			b = VecAppend(b, VecFromIdx(-lb, arglb));
		}
	}

	EVector argub;
	if (_ub != nullptr)
	{
		// upper bound
		const EVector& ub = *_ub;
		ASSERT(ub.rows() == numVars);
		argub = findUb(ub);
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

	EVector X;
	if (_x0 != nullptr) { X = *_x0; }
	else { X = EVector(numVars); X.setZero(); }

	EMatrix Z; // null space.
	
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

		EMatrix Q, R;
		EigenQrDecomp(Aeq.transpose(), &Q, &R, nullptr);
		if (numEqCstr < numVars)
			Z = MatFromColIdx(Q, colon(numEqCstr, numVars - 1)); 

		// is this necessary?
		{
			// find a feasible point for equality constraints
			EVector diff = Aeq * X - beq;
			if (VecAbs(diff).maxCoeff() > tolCons)
				X = EigenColPivQrSolve(Aeq, beq);
		}

		// is this necessary?
		if (numEqCstr > numVars)
		{
			EVector diff = Aeq * X - beq;
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
				EVector diff2 = A * X - b;
				if (diff2.maxCoeff() > 1e-8)
				{
					// exiting: the constraints or bounds are overly stringent
					// there's no feasible solution
					// equality constraints have been met
					exitFlag = -1;
				}
			}

			CalcInfeasibleLambda(isqp, lls, Q, R, H, f, X, eqix, indepIdx, normf, normA, lambda);
			return exitFlag;
		}

		if (Z.cols() == 0)
		{
			// Aeq's null space is null, there's no space for X to move.
			// and also X is infeasible.
			exitFlag = res.exitFlag;
			CalcInfeasibleLambda(isqp, lls, Q, R, H, f, X, eqix, indepIdx, normf, normA, lambda);
			EVector diff2 = A * X  - b;
			if (diff2.maxCoeff() > 1e-8)
			{
				// exiting: the constraints or bounds are overly stringent
				// there's no feasible solution
				// equality constraints have been met
				exitFlag = -1;
			}

			return exitFlag;
		}

		// check whether in Phase 1 of feasiblity point finding.
		// this part is removed.
		// if (verbosity == -2) ...
	}
	else
	{
		Z = EMatrix(1, 1); Z(0, 0) = 1;
	}

	// find initial feasible solution
	if (numCstr > 0 && numEqCstr < numCstr)
	{
		EVector cstr = A * X - b;
		float mc = VecFromIdx(cstr, colon(numEqCstr, numCstr - 1)).maxCoeff();
		if (mc > eps)
		{
			// add one slack variable and solve a linear system to get init feasible solution.
			int rowsA = A.rows();
			int colsA = A.cols();
			EMatrix A2(rowsA + 1, colsA + 1);
			A2.block(0, 0, rowsA, colsA) = A;
			A2.block(rowsA, 0, 1, colsA).setZero();
			A2.block(0, colsA, numEqCstr, 1).setZero();
			A2.block(numEqCstr, colsA, numCstr - numEqCstr + 1, 1).setOnes();
			A2.block(numEqCstr, colsA, numCstr - numEqCstr + 1, 1) *= -1.f;
			int quiet = -2;

			EMatrix feasM;
			EVector feasF(numVars); feasF.setZero();
			EVector vt(1); vt(0) = 1e-5;
			EVector b2 = VecAppend(b, vt);
			EVector vmc(1); vmc(0) = mc+ 1;
			EVector feasX = VecAppend(X, vmc);

			qpsub(feasM, feasF, // obj func
				A2, b2,  // constrs
				nullptr, nullptr, // lb & ub
				&feasX, numEqCstr, 
				kQpsub, eps);
		}
	}

	return 0;
}
