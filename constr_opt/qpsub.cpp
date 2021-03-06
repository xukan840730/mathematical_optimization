#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../common/maybe.h"
#include "qpsub.h"
#include "eqnsolv.h"
#include "compdir.h"

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

static EVector CalcLambda1(bool isqp, bool lls, const EMatrix& Q, const EMatrix& R, const EMatrix& H, const EVector& f, const EVector& X)
{
	EVector actLambda;
	if (isqp)
		actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H*X + f));
	else if (lls)
		actLambda = EigenColPivQrSolve(-R, Q.transpose()*(H.transpose() * (H*X - f)));
	else
		actLambda = EigenColPivQrSolve(-R, Q.transpose() * f);
	return actLambda;
}

static void CalcInfeasibleLambda(bool isqp, bool lls, const EMatrix& Q, const EMatrix& R, const EMatrix& H, const EVector& f, const EVector& X, const EVector& eqix, const EVector& indepIdx, 
	float normf, const EVector& normA, EVector& lambda)
{
	EVector actLambda = CalcLambda1(isqp, lls, Q, R, H, f, X);
	const EVector normedActLambda = normf * VecDivVec(actLambda, VecFromIdx(normA, eqix));
	VecReplaceRows(lambda, normedActLambda, VecFromIdx(indepIdx, eqix));
}

qpsubres qpsub(const EMatrix& H, const EVector& _f, const EMatrix& _A, const EVector& _b, const EVector* _lb, const EVector* _ub, const EVector* _x0, int numEqCstr, QpsubCaller caller, float eps) 
{
	ASSERT(_A.rows() == _b.rows());
	if (_x0 != nullptr) ASSERT(_x0->rows() == H.cols() || H.cols() == 0);
	int exitFlag = 1;
	int numIters = 0;

	// copy original matrix and vector.
	EVector f = _f;
	EMatrix A = _A;
	EVector b = _b;

	bool verbosity = false;

	ASSERT(numEqCstr <= A.rows());

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
	const float tolCons = 1e-10f;
	EVector lambda(numCstr); lambda.setZero();
	EVector aix = lambda;
	EVector eqix = colon(0, numEqCstr - 1);
	EVector indepIdx = colon(0, numCstr - 1); // independent constraint indices
	
	EMatrix Q(numVars, numVars); Q.setZero();
	EMatrix R;

	EVector X;
	if (_x0 != nullptr) { X = *_x0; }
	else { X = EVector(numVars); X.setZero(); }

	EMatrix Z; // null space.
	
	EVector actSet;
	EVector actInd;
	int actCnt = 0;	

	if (numEqCstr > 0)
	{
		eqnres res = eqnsolv(A, b, eqix, numVars, eps);
		if (res.rmvIdx.rows() > 0)
		{
			indepIdx = VecRmvIdx(indepIdx, res.rmvIdx);
			normA = VecFromIdx(normA, indepIdx);
		}
		actSet = MatFromRowIdx(A, eqix);
		actInd = eqix;
		actCnt = eqix.rows();

		if (res.exitFlag == -1)
		{
			// equalities are inconsistent, so x and lamabda have no valid values
			// return original x and zero for lambda
			return qpsubres(res.exitFlag, X, lambda, numIters);
		}

		// update flags.
		numEqCstr = eqix.rows();
		numCstr = A.rows();

		{
			EVector a(numEqCstr); a.setOnes();
			EVector b(numCstr - numEqCstr); b.setZero();
			aix = VecAppend(a, b);
		}
		const EMatrix Aeq = MatFromRowIdx(A, eqix);
		const EVector beq = VecFromIdx(b, eqix);

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
				return qpsubres(exitFlag, X, lambda, numIters);
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
			return qpsubres(exitFlag, X, lambda, numIters);
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

			return qpsubres(exitFlag, X, lambda, numIters);
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
	EVector cstr = A * X - b;
	if (numCstr > 0 && numEqCstr < numCstr)
	{
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
			EVector feasF(numVars + 1); feasF.setZero(); feasF(numVars) = 1;
			EVector vt(1); vt(0) = 1e-5f;
			EVector b2 = VecAppend(b, vt);
			EVector vmc(1); vmc(0) = mc+ 1;
			EVector feasX = VecAppend(X, vmc);

			qpsubres res2 = qpsub(feasM, feasF, // obj func
				A2, b2,  // constrs
				nullptr, nullptr, // lb & ub
				&feasX, numEqCstr, 
				kQpsub, eps);

			X = res2.X.block(0, 0, numVars, 1);
			cstr = A * X - b;
			float slackVar = res2.X(numVars);
			if (slackVar > eps)
			{
				exitFlag = -1;
				// exiting: the constraints are olvery stringent
				//lambda 
				const EVector lambdaS = res2.lambda.block(0, 0, numCstr, 1);
				const EVector t2 = normf * VecDivVec(lambdaS, normA);
				VecReplaceRows(lambda, t2, indepIdx);
				return qpsubres(exitFlag, X, lambda, numIters);
			}
		}
	}

	EVector SD;
	SearchDir dirType;
	EVector gf;
	if (isqp)
	{
		printf("isqp\n");
		gf = H * X + f;	
		// SD=-Z*((Z'*H*Z)\(Z'*gf));
		CompDirRes res2 = CompDir(&Z, &H, &gf, numVars, &f, eps);
		SD = res2.SD;
		dirType = res2.dirType;
	}
	else if (lls)
	{
		printf("lls\n");
		ASSERT(false); // TODO:
	}
	else // lp
	{
		gf = f;
		SD = -Z * Z.transpose() * gf;
		dirType = SearchDir::kStpDesc;
		if (SD.norm() < 1e-10 && numEqCstr > 0)
		{
			// this happens when equality constraint is perpendicular 
			// to objective function
			EVector actLambda = EigenColPivQrSolve(-R, Q.transpose() * gf);
			const EVector t2 = normf * VecDivVec(actLambda, VecFromIdx(normA, eqix));
			const EVector rowIdx = VecFromIdx(indepIdx, eqix);
			VecReplaceRows(lambda, t2, rowIdx);
			return qpsubres(exitFlag, X, lambda, numIters);
		}
	}


	// the maximum number of iterations for a simplex type method is when ncstr >= n:
	// maxiters = prod(1:ncstr)/(prod(1:numberOfVariables)*prod(1:max(1,ncstr-numberOfVariables)));

	// ------------------ main routine ----------------------------
	while (numIters <= maxIter)
	{
		numIters++;
		
		// find distance we can move in search direction SD before a constraint is violated.
		// gradient with respect to search direction.
		EVector GSD = A * SD;
		
		// note: we consider only constraints whose gradients are greater
		// than some threshold. If we considered all gradients greater than 
		// zero then it might be possible to add a constraint which would lead to
		// a singular (rank deficient) working set. The gradient (GSD) of such
		// a constraint in the direction of search would be very close to zero.

		EVector indf;
		EVector dist, ind2;
		int ind;
		float stepmin;
		{
			const EVector t2 = VecNot(aix);
			const EVector t3 = VecGt(GSD, errNorm * SD.norm());
			const EVector t4 = VecAnd(t3, t2);
			indf = findNnzRows(t4, 0);
		}
		if (indf.rows() == 0)
		{
			// no constraints to hit
			stepmin = 1e16f;
		}
		else
		{
			// find distance to the nearest constraint
			const EVector t2 = VecFromIdx(cstr, indf);
			const EVector t3 = VecFromIdx(GSD, indf);
			dist = VecAbs(VecDivVec(t2, t3));
			{
				stepmin = VecMin2(dist, ind2);
				float minRowIdx = VecMin(ind2);
				ind = indf(minRowIdx);
			}
		}
		// ---------- Update X ----------------
		
		// assume we do not delete a constraint
		bool delCstr = false;

		if (indf.rows() > 0 && ind2.rows() > 0)	// hit a constraint
		{
			if (dirType == SearchDir::kNewton) // newton step and hit a constraint
			{
				if (stepmin > 1) // overstepped minimum. reset stepmin
				{
					stepmin = 1;
					delCstr = true;
				}
				X = X + stepmin * SD;
			}
			else
			{
				// not a newton step and hit a constraint: is_qp or LLS or maybe lp
				X = X + stepmin * SD;
			}
		}
		else
		{
			// did not hit a constraint
			if (dirType == SearchDir::kNewton)  // newton step and no constraint hit: lls or maybe is_qp 
			{
				stepmin = 1; // exact distance to the solution. now delete constr.
				X = X + SD;
				delCstr = 1;
			}
			else // not a newton step: is_qp or lp or LLS
			{
				float smallRealEig = 0.f;
				if (isqp)
				{
					EMatrix ZHZ = Z.transpose() * H * Z;
					if (numVars < 400)
					{
						// if this is qp, i can reuse the small eigen value from compdir.
						smallRealEig = FindMinEigenVal(ZHZ); 
					}
					else
					{
						ASSERT(false);
					}
				}
				else
				{
					smallRealEig = 0;					
				}

				if ((!isqp && !lls) || (smallRealEig < -100*eps)) // LP or neg def: not LLS
				{
					// neg def -- unbounded
					if (SD.norm() > errNorm)
					{
						if (!normalize)
							stepmin = abs((X(numVars - 1) + 1e-5f) / (SD(numVars - 1)+eps));
						else
							stepmin = 1e16f;
						X = X + stepmin * SD;
						// how = 'unbounded'
						exitFlag = -1;
					}
					else // norm(SD) <= errNorm
					{
						// how = 'ill posed'
						exitFlag = -1;
					}

					if (verbosity)
					{
						if (SD.norm() > errNorm)
						{
							printf("Exiting: The solution is unbounded and at infinity;\n");
							printf("	 the constraints are not restrictive enough.\n");
						}
						else
						{
							printf("Exiting: The Search Direction is close to zero; \n");
							printf("	 the problem is ill-posed\n");
							printf("	 The gradient of the objective function may be zero\n");
							printf("	    or the problem may be badly conditioned.\n");
						}
					}
					return qpsubres(exitFlag, X, lambda, numIters);
				}
				else // singular: solve compatible system for a solution: isqp or lls
				{
					EMatrix projH;
					EMatrix projSD;
					EMatrix Zgf;
					if (isqp)
					{
						projH = Z.transpose() * H * Z;
						Zgf = Z.transpose() * gf;
						projSD = PseInv(projH) * (-Zgf);
					}
					else // LLS
					{ 
						ASSERT(false); // TODO:
					}

					// check if compatible
					if ((projH*projSD+Zgf).norm() > 10*eps*(projH.norm() + Zgf.norm()))
					{
						// system is incompatible --> it's a "chute": use SD from compdir
						// unbounded in SD direction
						if (SD.norm() > errNorm)
						{
							if (!normalize)
								stepmin = abs((X(numVars - 1) + 1e-5f) / (SD(numVars - 1) + eps));
							else
								stepmin = 1e16f;
							X = X + stepmin * SD;
							// how = 'unbounded'
							exitFlag = -1;
						}
						else // norm(SD) <= errNorm
						{
							// how = 'ill posed'
							exitFlag = -1;
						}

						if (verbosity)
						{
							if (SD.norm() > errNorm)
							{
								printf("Exiting: The solution is unbounded and at infinity;\n");
								printf("	 the constraints are not restrictive enough.\n");
							}
							else
							{
								printf("Exiting: The Search Direction is close to zero; \n");
								printf("	 the problem is ill-posed\n");
								printf("	 The gradient of the objective function may be zero\n");
								printf("	    or the problem may be badly conditioned.\n");
							}
						}
						return qpsubres(exitFlag, X, lambda, numIters);
					}
					else // convex -- move to minimum (compatible system)
					{
						SD = Z * projSD;
						dirType = SearchDir::kSingular;
						// first check if constraint is violated
						GSD = A * SD;
						{
							const EVector t2 = VecNot(aix);
							const EVector t3 = VecGt(GSD, errNorm * SD.norm());
							const EVector t4 = VecAnd(t3, t2);
							indf = findNnzRows(t4, 0);
						}

						if (indf.rows() == 0)
						{
							// no constraints to hit
							stepmin = 1;
							delCstr = true;
						}
						else
						{
							// find distance to the nearest constraint
							const EVector t2 = VecFromIdx(cstr, indf);
							const EVector t3 = VecFromIdx(GSD, indf);
							dist = VecAbs(VecDivVec(t2, t3));
							{
								stepmin = VecMin2(dist, ind2);
								float minRowIdx = VecMin(ind2);
								ind = indf(minRowIdx);
							}
						}
						if (stepmin > 1) // overstepped minimu; reset stepmin
						{
							stepmin = 1;
							delCstr = true;
						}
						X = X + stepmin * SD;
					}
				} // if ((!isqp && !lls) || (smallRealEig < -100*eps))
			} // if (dirType == SearchDir::kNewton)
		} // if (indf.rows() > 0 && ind2.rows() > 0)

		if (delCstr)
		{
			// note: only reach here if a minimum in the current subspace found
			if (actCnt > 0)
			{ 
				if (actCnt >= numVars - 1)
				{
					// avoid case when cind is greater than actCnt
					//if (cind <= actCnt)
					//{
					//	actSet = VecRmvIdx(actSet, cind);
					//	actInd = VecRmvIdx(actInd, cind);
					//}
				}
				EVector indLam;
				EVector rlambda;
				if (isqp)
				{
					rlambda = CalcLambda1(isqp, lls, Q, R, H, f, X);
					EVector actLambda = rlambda;
					VecReplaceRows(actLambda, VecAbs(VecFromIdx(rlambda, eqix)), eqix);
					
					indLam = findLtRows(actLambda, 0);
				}
				if (indLam.rows() == 0)
				{
					//lambda(indepInd(ACTIND)) = normf * (rlambda. / normA(ACTIND));
					const EVector normedRLambda = normf * VecDivVec(rlambda, VecFromIdx(normA, actInd));
					VecReplaceRows(lambda, normedRLambda, VecFromIdx(indepIdx, actInd));
					return qpsubres(exitFlag, X, lambda, numIters);
				}
				// remove constraint
				int lindn = VecMin3(actInd, indLam);
				EVector lind(1); lind(0) = lindn;
				VecRmvIdx(actSet, lind);
				aix(actInd(lindn)) = 0;
				EigenQrDelete(Q, R, lind);
				VecRmvIdx(actInd, lind);
				actCnt -= 2;
				simplexIter = false;
				ind = 0;
			}
			else // actCnt == 0
			{
				return qpsubres(exitFlag, X, lambda, numIters);
			}

			delCstr = false;
		} // end of if (delCstr)

		// calculate gradient w.r.t objective at this point
		if (isqp)
		{
			gf = H * X + f;
		}
		else if (lls)
		{
			gf = H.transpose()*(H * X - f);
			// else gf=f still true
		}

	} // end of while (numIters <= maxIter)

	return qpsubres(0, X, lambda, numIters);
}
