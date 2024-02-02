#ifndef CUSER_DATA_PAR_H
#define CUSER_DATA_PAR_H

#include "CUserData.h"

class CLevelPar;

typedef REAL* (*OLD_SOL_PAR)(const REAL_D);
typedef void (*BUILD_RHS_TO_PASS_PAR)(const EL_INFO*, const FE_SPACE**, REAL**,
		OLD_SOL_PAR, REAL, REAL);
typedef void (*BUILD_STATIC_MATRIX_TO_PASS_PAR)(const FE_SPACE**, DOF_MATRIX***, OLD_SOL_PAR, REAL, REAL);
typedef void (*BUILD_DYNAMIC_MATRIX_TO_PASS_PAR)(const FE_SPACE**, DOF_MATRIX***,
		OLD_SOL_PAR, REAL, REAL);
typedef REAL (*EXACT_SOL_PAR)(const REAL_D, REAL, int);
typedef REAL (*INIT_COND_PAR)(const REAL_D, REAL, int);
typedef REAL (*DIR_BOUND_COND_PAR)(const REAL_D, REAL, int);
typedef void (*MICRO_SOL_POSTPROCESS_PAR)(const FE_SPACE**, DOF_REAL_VEC**,
		REAL, REAL);

class CUserDataPar : public CUserData {
public:
	BUILD_STATIC_MATRIX_TO_PASS_PAR pBuildStatic; ///< pointer to the function which builds the dof matrix
	BUILD_DYNAMIC_MATRIX_TO_PASS_PAR pBuildDynamic; ///< pointer to the function which builds the dof matrix
	BUILD_RHS_TO_PASS_PAR pSetUpRHS; ///< pointer to the function which builds the right hand side of the system
	EXACT_SOL_PAR pExactSolution; ///< pointer to the function which returns exact solution (if testing)
	INIT_COND_PAR pInitialCondition; ///< pointer to the function which returns exact solution (if testing)
	DIR_BOUND_COND_PAR pDirBoundCond; ///< pointer to the function which returns the prescribed dirichlet data
	MICRO_SOL_POSTPROCESS_PAR pMicroSolPostprocess;
	int IsBiLinFormDependentOnPrevStep;
	double StartTime;
	double EndTime;
	double TimeStep;
	int NumTimeStep;
	void LoadUserDataFromFile(const char * pFileName);
	int CheckUserData();
	CUserDataPar() :
		pBuildStatic(NULL), pBuildDynamic(NULL), pSetUpRHS(NULL), pExactSolution(NULL),
				pInitialCondition(NULL), pDirBoundCond(NULL),
				pMicroSolPostprocess(NULL) {
	}
};

#endif //#ifndef C_USER_DATA_PAR_H
