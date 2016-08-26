#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "compdir.h"
#include <Eigen/Dense>

CompDirRes CompDir(const EMatrix* _Z, const EMatrix* _H, const EVector* _gf, int numVars, const EVector* _f, float eps)
{
	// CompDir computes a search direction in a subspace defined by Z.
	// return Newton direction if possible
	// return Random direction if gradient is small
	// otherwise, returns steepest descent direction
	// If the steepest descent direction is small it computes a negative
	// curvature direction based on the most negative eigenvalue.
	// For singular matrices, returns steepest descent even if small.
	const EMatrix& Z = *_Z;
	const EMatrix& H = *_H;
	const EVector& gf = *_gf;
	const EVector& f= *_f;

	EVector SD;
	SearchDir dirType;
	// SD=-Z*((Z'*H*Z)\(Z'*gf));
	// compute the projected newton direction if possible
	EMatrix projH = Z.transpose() * H * Z;

	EMatrix R;
	bool valid = EigenLlt(projH, &R);
	if (valid)
	{
		// positive definite: use Newton direction
		EVector t1 = Z.transpose() * gf;
		EVector t2 = EigenColPivQrSolve(R.transpose(),  t1);
		EVector t3 = EigenColPivQrSolve(R, t2);
		SD = -Z * t3;
		dirType = SearchDir::kNewton;
	}
	else
	{
		// not positive definite
		// if the gradient is small, try a random direction:
		// sometimes the search direction goes to zero in negative definite problems 
		// when current point rests on the top of quadratic function. 
		// In this case we can move in any direction to get an improvement in the function 
		// so that foil search by giving a random gradient.
		if (gf.norm() < sqrt(eps))
		{
			SD = -Z * Z.transpose() * MatAddScalar(MatRand(numVars, 1), - 0.5f);
			dirType = SearchDir::kRandom;
		}
		else 
		{
			// steepest descent
			const EVector stpDesc = -Z * (Z.transpose() * gf);
			// check if ||SD|| is close to zero
			if (stpDesc.norm() > sqrt(eps))
			{
				SD = stpDesc;
				dirType = SearchDir::kStpDesc;
			}
			else
			{
				// look for a negative curvature direction
				// some attempt at efficiency: usually it's
				// faster to use EIG unless many variables.
				EVector ev;
				if (numVars < 400)
				{
					// Matrix P is unbound in null-space of constraint matrix A.
					EMatrix::EigenvaluesReturnType eigvalP;
					Eigen::EigenSolver<EMatrix>::EigenvectorsType eigvecP;
					EigenValVec(projH, &eigvalP, &eigvecP);

					float minEigVal;
					int indminR = FindMinEigenValIdx(&eigvalP, &minEigVal);

					if (minEigVal < 0)
					{
						// check the sign of SD and the magnitude
						float SDtol = 100 * eps* gf.norm(); // we know gf.norm() > sqrt(eps)
						EVector Zev = Z * ev;
						float dotProd = Zev.dot(gf);
						if (dotProd > SDtol)
						{
							SD = -Zev;
							dirType = SearchDir::kEigenvec;
						}
						else if (dotProd < SDtol)
						{
							SD = Zev;
							dirType = SearchDir::kEigenvec;
						}
						else
						{
							SD = stpDesc;
							dirType = SearchDir::kStpDesc;
						}
					}
					else
					{
						// the projected hessian is singular
						// will propagate through the algorithm
						SD = stpDesc;
						dirType = SearchDir::kStpDesc;
					}
				}
				else
				{
					 //TODO: for complex or unsymmetric problem only.
					 ASSERT(false);
				}
			}
		}
	}

	CompDirRes res;
	res.SD = SD;
	res.dirType = dirType;
	return res;
}

