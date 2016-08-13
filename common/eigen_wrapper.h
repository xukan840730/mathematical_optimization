// a file to place eigen library function, it's huge and slow.

#ifndef _EIGEN_WRAPPER_H_
#define _EIGEN_WRAPPER_H_

void ChangeEVector(EVector* inout, int numRows);
void ChangeEMatrix(EMatrix* inout, int numRows);

EMatrix MatrixFromRowIdx(const EMatrix& A, const EVector& rowArr);
EMatrix MatrixFromColIdx(const EMatrix& A, const EVector& colArr);

// get row indices array after removing rows 
EVector RmRows(int numRows, const EVector& rmRow);
EMatrix MatrixRmvRowIdx(const EMatrix& A, const EVector& rmRow);
EVector VectorFromIdx(const EVector& a, const EVector& rowArr);
EVector VectorRmvIdx(const EVector& A, const EVector& rmIdx);

bool anyNonzero(const EMatrix& A, float eps);
// find non zeros and store them in row and col index.
void findNonzeros(const EMatrix& m, EVector* rowIdx, EVector* colIdx);
EVector findZeroIdx(const EVector& a, float eps);
EVector findNnzIdx(const EVector& a, float eps);
EVector colon(int j, int k); // matlab : operator
EVector VecAppend(const EVector& a, const EVector& b);

// eigen operation
bool EigenLlt(const EMatrix& m, EMatrix* l);

// TODO: to avoid brining Eigenvalue and Eigenvector types here.
//void EigenValVec(const EMatrix& m, 
//	EMatrix::EigenvaluesReturnType* eigenval, 
//	Eigen::EigenSolver<EMatrix>::EigenvectorsType* eigenvec);

void EigenValVec(const EMatrix& m, 
	void* eigenval, 
	void* eigenvec);

void EigenQrDecomp(const EMatrix& m, EMatrix* q, EMatrix* r, EMatrix* p = nullptr);

EVector EigenColPivQrSolve(const EMatrix& A, const EVector& b);

#endif
