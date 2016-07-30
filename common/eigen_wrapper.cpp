#include "common_shared.h"
#include "eigen_wrapper.h"

//-----------------------------------------------------------------------------------------------------------//
// helper function for eigen vectors.
//-----------------------------------------------------------------------------------------------------------//
void ChangeEVector(EVector* inout, int numRows)
{
	int numORows = inout->rows();
	if (numORows < numRows)
	{
		inout->conservativeResize(numRows);
		(*inout).block(numORows, 0, numRows - numORows, 1).setZero();
	}
	else if (numORows > numRows)
	{
		inout->conservativeResize(numRows);
	}
}

void ChangeEMatrix(EMatrix* inout, int numRows)
{
	int numORows = inout->rows();
	if (numORows < numRows)
	{
		inout->conservativeResize(numRows, numRows);

		(*inout).block(numORows, 0, numRows - numORows, numRows - numORows).setZero();
		(*inout).block(0, numORows, numRows - numORows, numRows - numORows).setZero();
		(*inout).block(numORows, numORows, numRows - numORows, numRows - numORows).setZero();
	}
	else
	{
		inout->conservativeResize(numRows, numRows);
	}
}

//------------------------------------------------------------------------------------------------------------//
bool EigenLlt(const EMatrix& m, EMatrix* l)
{
	Eigen::LLT<EMatrix> lltOfM; 
	lltOfM.compute(m);
	bool isPosD = lltOfM.info() == Eigen::Success;	// whether H is positive definite matrix.
	if (isPosD)
		*l = lltOfM.matrixL();

	return isPosD;
}


void EigenValVec(const EMatrix& m, 
	EMatrix::EigenvaluesReturnType* eigenval, 
	Eigen::EigenSolver<EMatrix>::EigenvectorsType* eigenvec)
{
	Eigen::EigenSolver<EMatrix> eh(m);
	*eigenval = eh.eigenvalues();
	*eigenvec = eh.eigenvectors();
}


void EigenQrDecomp(const EMatrix& m, EMatrix* q, EMatrix* r)
{
	Eigen::HouseholderQR<EMatrix> qrOfM(m);
	*q = qrOfM.householderQ();
	*r = qrOfM.matrixQR().triangularView<Eigen::Upper>();
}

EVector EigenColPivQrSolve(const EMatrix& A, const EVector& b)
{
	// solve A.x = b
	ASSERT(A.rows() == b.rows());
	EVector x = A.colPivHouseholderQr().solve(b);
	return x;
}
