//-----------------------------------------------------------------------------------------------//
// define common stuffs will be shared by all projects

#ifndef _COMMON_SHARED_H_
#define _COMMON_SHARED_H_

#define xassert(f) { if (!(f)) { int* p = 0; *p = 0; } }
#define ASSERT(f) { if (!(f)) { int* p = 0; *p = 0; } }

#define NDI_FLT_EPSILON		(1.19209290e-07F)
#define NDI_FLT_MAX			(3.40282347e+38f)
#define NDI_FLT_MIN			(1.17549435e-38F)

#include "memory.h"
#include "math.h"

#include <Eigen/Dense>

#define ND_MAX_EIGENSIZE 32
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, 0, 24, ND_MAX_EIGENSIZE> EMatrix;	// maximum 24 row, 32 cols
typedef Eigen::Matrix<float, Eigen::Dynamic, 1, 0, ND_MAX_EIGENSIZE, 1> EVector;

//-----------------------------------------------------------------------------------------------//
// function definition: 
// @input an vector of parameters
// @return scalar value.
//-----------------------------------------------------------------------------------------------//
typedef std::function<float(const EVector& input)> ScalarFunc;

//-----------------------------------------------------------------------------------------------//
// gradient function definition:
// @input a vector of parameters
// @output the gradient vector.
//-----------------------------------------------------------------------------------------------//
typedef std::function<void(const EVector& input, EVector* output)> GradientFunc;

//-----------------------------------------------------------------------------------------------//
// hessian function definition:
// @input a vector of parameters
// @output the hessian matrix
//-----------------------------------------------------------------------------------------------//
typedef std::function<void(const EVector& input, EMatrix* output)> HessianFunc;

//-----------------------------------------------------------------------------------------------//
// Continuous Differentiable function
// For 1 time differentiable function, we can get its gradient function.
// For 2 times differentiable function, we can get its hessian function
//-----------------------------------------------------------------------------------------------//
struct CD1Func {
	const ScalarFunc& f;
	const GradientFunc& g;

	CD1Func(const ScalarFunc& _f, const GradientFunc& _g) : f(_f) , g(_g) {}
};

struct CD2Func {
	const ScalarFunc& f;
	const GradientFunc& g;
	const HessianFunc& h;

	CD2Func(const ScalarFunc& _f, const GradientFunc& _g, const HessianFunc& _h) : f(_f) , g(_g) , h(_h) {}
};


#endif
