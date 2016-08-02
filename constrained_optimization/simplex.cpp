#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../common/bit_array.h"

typedef Eigen::Matrix<int, Eigen::Dynamic, 1, 0, ND_MAX_EIGENSIZE, 1> BasicVarIdx;

void MatrixRowDiv(EMatrix& m, int rowIdx, float f)
{
	m.row(rowIdx) /= f;
}

void Simplex1(EMatrix& tableu, const BasicVarIdx& basicVarIdx)
{
	ASSERT(tableu.cols() == basicVarIdx.rows());
	int numConstrs = tableu.rows() - 1;
	int numVars = tableu.cols() - 1;

	// scale the basic variables to 1
	for (int ii = 0; ii < numVars; ii++)
	{
		int colIdx = basicVarIdx(ii);
		if (colIdx >= 0)
		{
			ASSERT(colIdx < numConstrs);
			float scale = tableu(colIdx, ii);
			ASSERT(scale != 0.f);

			MatrixRowDiv(tableu, colIdx, scale);
		}
	}
}


void Simplex(EMatrix& tableu, BasicVarIdx& basicVarIdx)
{
	// find a non basic variable with negative lambda.
	ASSERT(tableu.cols() == basicVarIdx.rows());
	int numConstrs = tableu.rows() - 1;
	int numVars = tableu.cols() - 1;

	float minLambda = 0.f;
	int minLambdaIdx = -1;
	for (int ii = 0; ii < numVars; ii++)
	{
		int colIdx = basicVarIdx(ii);
		if (colIdx < 0) // non basic variable
		{
			float lambda = tableu(numConstrs, ii);
			if (lambda < minLambda)
			{
				minLambda = lambda;
				minLambdaIdx = ii;
			}
		}
	}

	if (minLambdaIdx < 0) 
	{
		// algorithm termiate.
	}
	else
	{
		// find a row with minimum ratio	
	}
}
