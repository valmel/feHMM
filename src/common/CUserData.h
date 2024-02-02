#ifndef C_USER_DATA_H
#define C_USER_DATA_H
#include <iostream>
/// CLevelData is structure holding all info about the problem supported by user

#ifndef _ALBERTA_H_
#include "alberta.h"
#endif

#ifdef USING_MPI
#include "mpi.h"
#endif

#include "ErrorReturnCodes.h"

typedef void (*BUILD_RHS_TO_PASS)(const EL_INFO*, const FE_SPACE**, REAL**);
typedef void (*BUILD_STATIC_MATRIX_TO_PASS)(const FE_SPACE**, DOF_MATRIX***);
typedef void (*BUILD_DYNAMIC_MATRIX_TO_PASS)(const FE_SPACE**, DOF_MATRIX***);
typedef void (*MICRO_SOL_POSTPROCESS)(const FE_SPACE**, DOF_REAL_VEC**);
typedef REAL (*EXACT_SOL)(const REAL_D, int);
typedef REAL (*DIR_BOUND_COND)(const REAL_D, int);

struct SolverInfo {
	REAL tol;
	int miter;
	int info;
	int icon;
	int restart;
	OEM_SOLVER solver;
};

class CUserData {
protected:
	void ExtractParamFromLine(const char*, int, int, void *); ///< function splits line in case of several parameters
public:
#ifdef USING_MPI
	int Rank;
#endif
	int ProblemNum; ///< number of problem user wants to solve
	int NumOfLevels; ///< number of levels, i.e. number of scales
	int DimOfWorld; ///< dimension of the world, i.e. if f(x,y) then DimOfWorld=2
	int DimOfProb; ///< dimension of the problem, i.e. if f=(f1,f2,f3) then DimOfProb=3
	int *pNumOfRefin; ///< numbers of refinements for each level
	int *pQuadDeg; ///< degree of quadrature for the lowest level (but u can decide later there ;))
	REAL *pScaleFactor; ///< vector of relative dimensions of microcubes to the previous mesh element dimension (only regular meshes supported)
	int *pFEdeg; ///< degree of lagrange finite elements on the lowest level
	int **pMatrixMask; ///< at each level says which parts of ProblemMatrix should ne aasembled - DimOfProb x DimOfProb
	BUILD_STATIC_MATRIX_TO_PASS pBuildStatic; ///< pointer to the function which builds the dof matrix
	BUILD_DYNAMIC_MATRIX_TO_PASS pBuildDynamic; ///< pointer to the function which builds the dof matrix
	BUILD_RHS_TO_PASS pSetUpRHS; ///< pointer to the function which builds the right hand side of the system
	EXACT_SOL pExactSolution; ///< pointer to the function which returns exact solution (if testing)
	DIR_BOUND_COND pDirBoundCond; ///< pointer to the function which returns the prescribed dirichlet data
	MICRO_SOL_POSTPROCESS pMicroSolPostprocess; ///< pointer to the function which postprocesses solution on micro level
	char *pMacroFile; ///< name of macro file obtaining the geometry (standard ALBERTA macro file)
	SolverInfo Solver;
	int BoundCond;
	CUserData() :
		NumOfLevels(-1), DimOfWorld(-1), DimOfProb(-1), pNumOfRefin(NULL),
				pQuadDeg(NULL), pScaleFactor(NULL), pFEdeg(NULL),
				pMatrixMask(NULL), pBuildStatic(NULL), pBuildDynamic(NULL),
				pSetUpRHS(NULL), pExactSolution(NULL), pDirBoundCond(NULL),
				pMicroSolPostprocess(NULL), pMacroFile(NULL), BoundCond(-1) {
	}///< data are "invalid" after constructor is called
	void LoadUserDataFromFile(const char *); ///< function extracts parameters from user provided init file
	int CheckUserData(); ///< testovanie uzivatelskych dat
	int CheckUserDataWithoutFuncPointers();
	void FillSolverData(SolverInfo &);
	~CUserData();
};
#endif // #ifndef C_USER_DATA_H
