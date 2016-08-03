#include "simplex.h"

void Simplex1(EMatrix& tableu, const BasicVarIdx& basicVarIdx)
{
	ASSERT(tableu.cols() == basicVarIdx.rows() + 1);
	int numConstrs = tableu.rows() - 1;
	int numVars = tableu.cols() - 1;

	for (int ii = 0; ii < numVars; ii++)
	{
		int rowIdx = basicVarIdx(ii);
		if (rowIdx >= 0)
		{
			// scale the basic variables to 1
			ASSERT(rowIdx < numConstrs);
			{
				float scale = tableu(rowIdx, ii);
				ASSERT(scale != 0.f);
				tableu.row(rowIdx) /= scale;
			}

			// elimate the lambda of a basic variable
			{
				float scale = tableu(numConstrs, ii);
				if (scale != 0.f)
					tableu.row(numConstrs) -= tableu.row(rowIdx) * scale;
			}
		}
	}
}


int Simplex(EMatrix& tableu, BasicVarIdx& basicVarIdx)
{
	// find a non basic variable with negative lambda.
	ASSERT(tableu.cols() == basicVarIdx.rows() + 1);
	int numConstrs = tableu.rows() - 1;
	int numVars = tableu.cols() - 1;

	Simplex1(tableu, basicVarIdx);

	int res = 0;

	while (true)
	{
		float minLambda = 0.f;
		int newBasicIdx = -1;	// a non-basic variable will become basic.
		for (int ii = 0; ii < numVars; ii++)
		{
			int colIdx = basicVarIdx(ii);
			if (colIdx < 0) // this is non basic variable
			{
				float lambda = tableu(numConstrs, ii);
				if (lambda < minLambda)
				{
					minLambda = lambda;
					newBasicIdx = ii;
				}
			}
		}

		if (newBasicIdx < 0) 
		{
			// algorithm termiate.
			res = 0;
			break;
		}
		else
		{
			// find a row with minimum ratio
			float minRatio = NDI_FLT_MAX;
			int replacedRowIdx = -1;
			for (int ii = 0; ii < numConstrs; ii++)
			{
				float coeff = tableu(ii, newBasicIdx);
				if (coeff > 0.f)
				{
					float ratio = tableu(ii, numVars) / coeff;
					if (ratio < minRatio)
					{
						minRatio = ratio;
						replacedRowIdx = ii;
					}
				}
			}

			// swap basic and non-basic variables here.
			if (replacedRowIdx < 0)
			{
				// unbounded.
				res = 1;
				break;
			}
			else
			{
				{
					float scale = tableu(replacedRowIdx, newBasicIdx);
					tableu.row(replacedRowIdx) /= scale;
				}

				for (int ii = 0; ii < tableu.rows(); ii++)
				{
					if (ii != replacedRowIdx)
					{
						float scale = tableu(ii, newBasicIdx);
						tableu.row(ii) -= tableu.row(replacedRowIdx) * scale;
					}
				}

				// update basic variable flags
				for (int ii = 0; ii < numVars; ii++)
				{
					if (basicVarIdx(ii) == replacedRowIdx)
					{
						basicVarIdx(newBasicIdx) = basicVarIdx(ii);
						basicVarIdx(ii) = -1;
						break;
					}
				}
			}
		}
	}

	printf("Simplex done! %d\n", res);
	return res;
}

int SolveInitFeasible(const EMatrix& A, const EVector& b, EVector* x0)
{
	ASSERT(A.rows() == b.rows());
	ASSERT(x0->rows() == A.cols());
	int numVars = A.cols();
	int numConstrs = A.rows();

	int numRows = numConstrs + 1;
	int numCols = numVars * 2 + numConstrs * 2 + 1;

	EMatrix tableu(numRows, numCols);

	for (int ii = 0; ii < numConstrs; ii++)
	{
		for (int jj = 0; jj < numVars; jj++)
		{
			float aij = A(ii, jj);
			tableu(ii, jj) = aij;
			tableu(ii, numVars + jj) = -aij;
		}
	}

	tableu.block(0, numVars * 2, numConstrs, numConstrs).setIdentity();
	tableu.block(0, numVars * 2 + numConstrs, numConstrs, numConstrs).setIdentity();
	tableu.block(0, numVars * 2 + numConstrs * 2, numConstrs, 1) = b;

	tableu.block(numConstrs, 0, 1, numVars * 2).setZero();
	tableu.block(numConstrs, numVars * 2, 1, numConstrs).setOnes();
	tableu.block(numConstrs, numVars * 2 + numConstrs, 1, numConstrs).setZero();
	tableu(numRows - 1, numCols - 1) = 0;

	BasicVarIdx bvarIdx(numVars * 2 + numConstrs * 2);
	for (int ii = 0; ii < bvarIdx.rows(); ii++)
		bvarIdx(ii) = -1;
	for (int ii = 0; ii < numConstrs; ii++)
		bvarIdx(numVars * 2 + ii) = ii;

	int res = Simplex(tableu, bvarIdx);
	if (res == 0)
	{
		float f = tableu(numRows - 1, numCols - 1);
		if (fabs(f) < NDI_FLT_EPSILON)	// if initial feasible solution is found, the sum of artifical variables must be zero
		{
			EVector initVars(numVars * 2 + numConstrs * 2);
			for (int ii = 0; ii < initVars.rows(); ii++)
			{
				initVars(ii) = bvarIdx(ii) >= 0 ? tableu(bvarIdx(ii), numCols - 1) : 0.f;
			}

			// convert from, x = y - z.
			for (int ii = 0; ii < numVars; ii++)
				(*x0)(ii) = initVars(ii) - initVars(ii + numVars);
		}
		else
		{
			res = 1;
		}
	}

	return res;
}
