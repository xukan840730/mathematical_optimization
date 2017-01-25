// a file to place eigen library function, it's huge and slow.

#ifndef _EIGEN_WRAPPER_H_
#define _EIGEN_WRAPPER_H_

void ChangeEVector(EVector* inout, int numRows);
void ChangeEMatrix(EMatrix* inout, int numRows);

EMatrix MatFromRowIdx(const EMatrix& A, const EVector& rowArr);
EMatrix MatFromColIdx(const EMatrix& A, const EVector& colArr);

// get row indices array after removing rows 
EVector RmRows(int numRows, const EVector& rmRow);
EVector RmCols(int numCols, const EVector& rmCol);
EMatrix MatRmvRowIdx(const EMatrix& A, const EVector& rmRow);
EMatrix MatRmvColIdx(const EMatrix& A, const EVector& rmCol);
EVector VecFromIdx(const EVector& a, const EVector& rowArr);
EVector VecRmvIdx(const EVector& A, const EVector& rmIdx);
void VecReplaceRows(EVector& a, const EVector& b, const EVector& rowIdx);

bool anyNnz(const EMatrix& A, float eps);
// find non zeros and store them in row and col index.
typedef bool (*condf)(float x, void* params);  // a condition function takes float as parameter
EVector findRows(const EVector& a, condf t, void* params);  // find all rows satisfy condition func

void findNonzeros(const EMatrix& m, EVector* rowIdx, EVector* colIdx);
EVector findZeroRows(const EVector& a, float eps);
EVector findNnzRows(const EVector& a, float eps);
EVector findLtRows(const EVector& a, float b);
EVector findGtRows(const EVector& a, float b);

EVector VecCond(const EVector& a, condf t, void* params);  // i don't know a better name for this func
EVector VecEq(const EVector& a, float f);
EVector VecGt(const EVector& a, float f);	// greater than
EVector VecSt(const EVector& a, float f);	// smaller than
EVector VecAbs(const EVector& a);
EVector VecDivVec(const EVector& a, const EVector& b);

// logic operation
EVector VecAnd(const EVector& a, const EVector& b);
EVector VecOr(const EVector& a, const EVector& b);
EVector VecNot(const EVector& a);

float VecMin(const EVector& a);
float VecMin2(const EVector& a, EVector& indices);
int VecMin3(const EVector& a, const EVector& rowIdx);

EVector colon(int j, int k); // matlab : operator
EVector VecAppend(const EVector& a, const EVector& b);
EMatrix MatRowAppend(const EMatrix& a, const EMatrix& b);
EMatrix MatColAppend(const EMatrix& a, const EMatrix& b);

EMatrix MatAddScalar(const EMatrix& m, float s);
EMatrix MatRand(int row, int col);

// eigen operation
bool EigenLlt(const EMatrix& m, EMatrix* l);

// TODO: to avoid brining Eigenvalue and Eigenvector types here.
//void EigenValVec(const EMatrix& m, 
//	EMatrix::EigenvaluesReturnType* eigenval, 
//	Eigen::EigenSolver<EMatrix>::EigenvectorsType* eigenvec);

void EigenValVec(const EMatrix& m, 
	void* eigenval, 
	void* eigenvec);
int FindMinEigenValIdx(void* eigenval, float* outMinEigenVal);
float FindMinEigenVal(const EMatrix& m);

void EigenQrDecomp(const EMatrix& m, EMatrix* q, EMatrix* r, EMatrix* p = nullptr);
void EigenQrDelete(EMatrix& q, EMatrix& r, const EVector& rmCols);

EVector EigenColPivQrSolve(const EMatrix& A, const EVector& b);

EMatrix PseInv(const EMatrix& e, float toler = 1e-6);

#endif
