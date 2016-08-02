#ifndef _AUGMENTED_LAGRANGIAN_H_
#define _AUGMENTED_LAGRANGIAN_H_

// Augmented Lagrangian method: this is in factor a penalty function like method, not Lagrangian method.
OptResult ALMethod(const CD1Func& objectiveF, const EVector& x0, int numEConstr, const CD1Func* econstrFs, const LagrangeMultMethodParams& params);
OptResult ALMethod(const CD1Func& objectiveF, const EVector& x0, int numEConstr, const CD1Func* econstrFs, int numInConstr, const CD1Func* inconstrFs, const LagrangeMultMethodParams& params);

#endif
