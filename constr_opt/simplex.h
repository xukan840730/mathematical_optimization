#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../common/bit_array.h"

typedef Eigen::Matrix<int, Eigen::Dynamic, 1, 0, ND_MAX_EIGENSIZE, 1> BasicVarIdx;

int Simplex(EMatrix& tableu, BasicVarIdx& basicVarIdx);
int SolveInitFeasible(const EMatrix& Aeq, const EVector& beq, const EMatrix& Ain, const EVector& bin, EVector* x0);

#endif
