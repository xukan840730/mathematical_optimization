#ifndef _QP_SUB_H_
#define _QP_SUB_H_

enum QpsubCaller
{
	kDefault,	// default
	kQpsub,		// qp sub problem
	kLsqlin,	// least square lin
};

void qpsub(const EMatrix& H, const EVector& _f, const EMatrix& A, const EVector& b, const EVector& lb, const EVector& ub, const EVector* x0, int numEqCstr, int numCstr, QpsubCaller caller); 

#endif
