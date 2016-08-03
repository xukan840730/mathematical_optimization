#include "simplex.h"

//-------------------------------------------------------------------------//
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

//-------------------------------------------------------------------------//
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

int SolveInitFeasible(const EMatrix* Aeq, const EVector* beq, const EMatrix* Ain, const EVector* bin, EVector* x0)
{
	if (Aeq || beq) ASSERT(Aeq && beq && Aeq->rows() == beq->rows());
	if (Ain || bin) ASSERT(Ain && bin && Ain->rows() == bin->rows());

	int numOVars = x0->rows(); // number of old variables.
	int numEconstrs = Aeq ? (*Aeq).rows() : 0;
	int numInconstrs = Ain ? (*Ain).rows() : 0;
	int numTotalConstrs = numEconstrs + numInconstrs;
	ASSERT(numTotalConstrs > 0);

	// each old varialbe is converted to xplus - xminus, and xplus >= 0, xminus >= 0,
	int numNRows = numEconstrs + numInconstrs + 1;
	int numNCols = numOVars * 2 + numEconstrs + numInconstrs * 2 + 1;

	EMatrix tableu(numNRows, numNCols);

	// convert Aeq.x to Aeq.(xplus - xminus)
	for (int ii = 0; ii < numEconstrs; ii++)
	{
		tableu.block(ii, 0, 1, numOVars) = (*Aeq).row(ii);
		tableu.block(ii, numOVars, 1, numOVars) = -(*Aeq).row(ii);
	}
	// convert Ain.x to Ain.(xplus - xminus)
	for (int ii = 0; ii < numInconstrs; ii++)
	{
		tableu.block(numEconstrs + ii, 0, 1, numOVars) = (*Ain).row(ii);
		tableu.block(numEconstrs + ii, numOVars, 1, numOVars) = -(*Ain).row(ii);
	}

	int teqColIdx = numOVars * 2;
	int sColIdx = teqColIdx + numEconstrs;
	int tinColIdx = sColIdx + numInconstrs;
	int bColIdx = tinColIdx + numInconstrs;
	ASSERT(bColIdx == numNCols - 1);

	// artifical variable teq, from Aeq.x + teq = beq
	tableu.block(0, teqColIdx, numEconstrs, numEconstrs).setIdentity();
	tableu.block(0, sColIdx, numInconstrs, numInconstrs * 2).setZero();

	// slack variable s, from Ain.x + s + tin = bin
	tableu.block(numEconstrs, teqColIdx, numEconstrs, numEconstrs).setZero();
	tableu.block(numEconstrs, sColIdx, numInconstrs, numInconstrs).setIdentity();
	//tableu.block(numEconstrs, sColIdx, numInconstrs, numInconstrs) *= -1.f;

	// artifical variable tin, from Ain.x + s + tin = bin
	tableu.block(numEconstrs, tinColIdx, numInconstrs, numInconstrs).setIdentity();

	// beq part
	if (beq != nullptr)
		tableu.block(0, bColIdx, numEconstrs, 1) = *beq;
	// bin part
	if (bin != nullptr)
		tableu.block(numEconstrs, bColIdx, numInconstrs, 1) = *bin;

	// fill the rest part
	tableu.block(numTotalConstrs, 0, 1, numOVars * 2).setZero();
	tableu.block(numTotalConstrs, teqColIdx, 1, numEconstrs).setOnes();
	tableu.block(numTotalConstrs, sColIdx, 1, numInconstrs).setZero();
	tableu.block(numTotalConstrs, tinColIdx, 1, numInconstrs).setOnes();
	tableu(numNRows - 1, numNCols - 1) = 0;

	int numNVars = numOVars * 2 + numEconstrs + numInconstrs * 2;
	BasicVarIdx bvarIdx(numNVars);
	for (int ii = 0; ii < numNVars; ii++)
		bvarIdx(ii) = -1;
	for (int ii = 0; ii < numEconstrs; ii++)
		bvarIdx(teqColIdx + ii) = ii;
	for (int ii = 0; ii < numInconstrs; ii++)
		bvarIdx(tinColIdx + ii) = numEconstrs + ii;

	int res = Simplex(tableu, bvarIdx);
	if (res == 0)
	{
		float f = tableu(numNRows - 1, numNCols - 1);
		if (fabs(f) < NDI_FLT_EPSILON)	// if initial feasible solution is found, the sum of artifical variables must be zero
		{
			EVector initVars(numNVars);
			for (int ii = 0; ii < numNVars; ii++)
			{
				initVars(ii) = bvarIdx(ii) >= 0 ? tableu(bvarIdx(ii), numNCols - 1) : 0.f;
			}

			// convert from, x = y - z.
			for (int ii = 0; ii < numOVars; ii++)
				(*x0)(ii) = initVars(ii) - initVars(ii + numOVars);
		}
		else
		{
			res = 1;
		}
	}

	return res;
}
