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

bool anyNnz(const EMatrix& A, float eps);
// find non zeros and store them in row and col index.
typedef bool (*condf)(float x, void* params);  // a condition function takes float as parameter
EVector findRows(const EVector& a, condf t, void* params);  // find all rows satisfy condition func

void findNonzeros(const EMatrix& m, EVector* rowIdx, EVector* colIdx);
EVector findZeroRows(const EVector& a, float eps);
EVector findNnzRows(const EVector& a, float eps);

EVector VecCond(const EVector& a, condf t, void* params);  // i don't know a better name for this func
EVector VecEq(const EVector& a, float f);

EVector colon(int j, int k); // matlab : operator
EVector VecAppend(const EVector& a, const EVector& b);
EMatrix MatRowAppend(const EMatrix& a, const EMatrix& b);
EMatrix MatColAppend(const EMatrix& a, const EMatrix& b);

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
