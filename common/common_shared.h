//-----------------------------------------------------------------------------------------------//
// define common stuffs will be shared by all projects

#ifndef _COMMON_SHARED_H_
#define _COMMON_SHARED_H_

#define xassert(f) { if (!(f)) { int* p = 0; *p = 0; } }

#define NDI_FLT_EPSILON		(1.19209290e-07F)
#define NDI_FLT_MAX			(3.40282347e+38f)
#define NDI_FLT_MIN			(1.17549435e-38F)

#include "memory.h"
#include "math.h"


#endif