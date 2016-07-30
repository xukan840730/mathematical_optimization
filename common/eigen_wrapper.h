// a file to place eigen library function, it's huge and slow.

#ifndef _EIGEN_WRAPPER_H_
#define _EIGEN_WRAPPER_H_

void ChangeEVector(EVector* inout, int numRows);
void ChangeEMatrix(EMatrix* inout, int numRows);


bool EigenLlt(const EMatrix& m, EMatrix* l);

void EigenValVec(const EMatrix& m, 
	EMatrix::EigenvaluesReturnType* eigenval, 
	Eigen::EigenSolver<EMatrix>::EigenvectorsType* eigenvec);

void EigenQrDecomp(const EMatrix& m, EMatrix* q, EMatrix* r);

EVector EigenColPivQrSolve(const EMatrix& A, const EVector& b);

#endif
