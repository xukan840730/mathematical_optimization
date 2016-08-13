#include "../common/common_shared.h"
#include "../common/eigen_wrapper.h"
#include "../linear_algebra/scalar_matrix.h"
#include "quad-prog.h"
#include <Eigen\Dense>

//------------------------------------------------------------------------------------------------------//
EQuadProgRes EQuadProg(const EMatrix& H, const EVector& q, const EMatrix& Aeq, const EVector& beq)
{
	ASSERT(Aeq.rows() > 0 && Aeq.cols() > 0);
	ASSERT(Aeq.rows() <= Aeq.cols());	// otherwise the linear system is overdetermined and could not have a solution.
	ASSERT(Aeq.rows() == beq.rows());
	ASSERT(H.rows() == H.cols());
	ASSERT(H.rows() == q.rows());
	ASSERT(Aeq.cols() == H.cols());

	EQuadProgRes result;
	result.type = EQuadProgRes::kUnbounded;

	EMatrix AT = Aeq.transpose();

	int mm = AT.rows();
	int nn = AT.cols();
	ASSERT(mm >= nn);

	// QR decomposition to get null space.
	EMatrix qOfAT;
	EMatrix rOfAT;
	EigenQrDecomp(AT, &qOfAT, &rOfAT, nullptr);

	EMatrix qHat = qOfAT.block(0, 0, mm, nn);
	EMatrix rHat = rOfAT.block(0, 0, nn, nn);

	// solve A * delta = B;
	EMatrix rHatT = rHat.transpose();
	EVector uu = EigenColPivQrSolve(rHatT, beq);
	bool slnEst = (rHatT * uu).isApprox(beq, 0.0001f);
	if (!slnEst)
		return result;

	EVector xHat = qHat * uu;
	if (mm == nn)
	{
		// the number of active constaints equals to the number of variables, the system has only 1 solution.
		result.type = EQuadProgRes::kFound;
		result.xstar = xHat;
		return result;
	}

	// let P = Qn^t . H . Qn, if P is positive definite, which means H is positive definite on the null space of constraint matrix A.
	// so we can get a unique solution.
	EMatrix qN = qOfAT.block(0, nn, mm, mm - nn);
	EMatrix P = qN.transpose() * H * qN;

	EMatrix L;
	bool isPosD = EigenLlt(P, &L);
	if (isPosD)
	{
		// TODO: replaced by Eigen library LLT.
		EMatrix PInv(L.rows(), L.cols());
		LLtInverse(&PInv, L);

		EMatrix ww = (xHat.transpose() * H * qN + q.transpose() * qN);
		EMatrix vv = -1.f * PInv * ww.transpose();

		EVector vPart = qN * vv;
		ASSERT(vPart.rows() == xHat.rows());
		EVector xstar = xHat + vPart;

		result.type = EQuadProgRes::kFound;
		result.xstar = xstar;
	}
	else
	{
		// Matrix P is unbound in null-space of constraint matrix A.
		EMatrix::EigenvaluesReturnType eigenvalP;
		Eigen::EigenSolver<EMatrix>::EigenvectorsType eigenvecP;
		EigenValVec(P, &eigenvalP, &eigenvecP);

		// find minimun eigenvec column.
		int indminR = 0;
		{
			ASSERT(nn == eigenvalP.rows());
			float smallestVal = eigenvalP(0).real();
			for (int ii = 1; ii < nn; ii++)
			{
				if (eigenvalP(ii).real() < smallestVal)
				{
					smallestVal = eigenvalP(ii).real();
					indminR = ii;
				}
			}
		}

		EVector eVrH = eigenvecP.col(indminR).real();
		ASSERT(eVrH.rows() == qN.cols());
		EVector stepp = qN * eVrH;

		result.type = EQuadProgRes::kUnbounded; 
		result.stepp = stepp;
	}

	return result;
}

