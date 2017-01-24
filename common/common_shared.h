//-----------------------------------------------------------------------------------------------//
// define common stuffs will be shared by all projects

#ifndef _COMMON_SHARED_H_
#define _COMMON_SHARED_H_

#define xassert(f) { if (!(f)) { int* p = 0; *p = 0; } }
#define ASSERT(f) { if (!(f)) { int* p = 0; *p = 0; } }

#define NDI_FLT_EPSILON		(1.19209290e-07F)
//#define NDI_FLT_MAX			(3.40282347e+38f)
#define NDI_FLT_MAX			(FLT_MAX)
#define NDI_FLT_MIN			(1.17549435e-38F)

#define ARRAY_COUNT(a) (sizeof(a) / sizeof(*(a)))

typedef long U32;
typedef long long U64;


static inline bool IsFinite(float a)
{
	return (*(unsigned int*)(&a) & 0x7F800000) != 0x7F800000;
}

#include "memory.h"
#include "math.h"

#include <Eigen/Core>

#define ND_MAX_EIGENSIZE 16 
typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, 0, ND_MAX_EIGENSIZE, ND_MAX_EIGENSIZE> EMatrix;	// maximum 24 row, 32 cols
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
public:
	const ScalarFunc* f;
	const GradientFunc* g;

public:
	CD1Func() : f(nullptr), g(nullptr) {}
	CD1Func(const ScalarFunc& _f, const GradientFunc& _g) : f(&_f), g(&_g) {}
};

struct CD2Func : public CD1Func {
	const HessianFunc* h;

public:
	CD2Func(const ScalarFunc& _f, const GradientFunc& _g, const HessianFunc& _h) : CD1Func(_f, _g), h(&_h) {}
};


struct LagrangeMultMethodParams
{
public:
	LagrangeMultMethodParams()
		: m_epsilon1(0.001f)
		, m_epsilon2(0.001f)
		, m_lamda1(1.f)
		, m_maxIter(20)
	{}

	float m_epsilon1;	// gradient epsilon to stop the iteration.
	float m_epsilon2;	
	float m_lamda1;		// init guess of lagrange multiplier
	int m_maxIter;		// max number of iterations
};

// optimization algorithm result.
struct OptResult
{
	EVector xstar;		// x which minimum 
	int numIter;		// number of iteration used.
};

//------------------------------------------------------------------------//
// Search Direction: defines which types of search direction
//------------------------------------------------------------------------//
enum SearchDir
{
	kNewton,
	kRandom,
	kEigenvec,
	kStpDesc,
};

#endif
