#ifndef _QP_SUB_H_
#define _QP_SUB_H_

enum QpsubCaller
{
	//kDefault,	// default
	kQpsub,		// qp sub problem
	kLsqlin,	// least square lin
};

struct qpsubres
{
	int exitFlag;
	EVector X;
	EVector lambda;

	qpsubres(int _exitFlag, const EVector& _X, const EVector& _lambda)
		: exitFlag(_exitFlag)
		, X(_X)
		, lambda(_lambda)
	{}
};

qpsubres qpsub(const EMatrix& H, const EVector& _f,
	const EMatrix& A, const EVector& b, 
	const EVector* lb, const EVector* ub, 
	const EVector* x0, int numEqCstr, 
	QpsubCaller caller, float eps); 


#endif
