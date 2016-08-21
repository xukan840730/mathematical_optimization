#ifndef _QP_SUB_H_
#define _QP_SUB_H_

enum QpsubCaller
{
	//kDefault,	// default
	kQpsub,		// qp sub problem
	kLsqlin,	// least square lin
};

int qpsub(const EMatrix& H, const EVector& _f,
	const EMatrix& A, const EVector& b, 
	const EVector* lb, const EVector* ub, 
	const EVector* x0, int numEqCstr, 
	QpsubCaller caller, float eps); 


#endif
