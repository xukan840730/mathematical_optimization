
#include "../linear_algebra/gradient.h"
#include "../linear_algebra/hessian.h"

void NewtonsMethod(const ScalarF F, const Gradient* g, const Hessian* H, const ScalarVector& initGuess, ScalarVector* result);