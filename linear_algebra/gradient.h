
#include "scalar_vector.h"

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

typedef std::function<void(const EVector& input, EVector* output)> GradientFunc;

#endif