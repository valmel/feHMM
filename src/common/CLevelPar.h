#ifndef CLEVEL_PAR_H
#define CLEVEL_PAR_H

#include "CLevel.h"
#include "CUserDataPar.h"
#include "CAlbertaInterfacePar.h"
#include "CTimePartitioning.h"
#include "SaveToVTKpar.h"

struct SUserInterfacePar {
	BUILD_STATIC_MATRIX_TO_PASS_PAR pBuildStatic;
	BUILD_DYNAMIC_MATRIX_TO_PASS_PAR pBuildDynamic;
	BUILD_RHS_TO_PASS_PAR pSetUpRHS;
	MICRO_SOL_POSTPROCESS_PAR pMicroSolPostprocess;
	EXACT_SOL_PAR pExact;
};

class CLevelPar: public CLevel{
	static REAL Time;
	static REAL UserLocalExact(const REAL_D);
	static EXACT_SOL_PAR pExact;
	REAL LocalExact(const REAL_D);

	static CTimePartitioning TimePar;
	int IsBiLinFormDependentOnPrevStep;
	SaveToVTKpar * pSaveVTKpar;
	static REAL** LocalDofVec;
	static REAL* OldSolution;
	static int DimOfProb;
#ifdef USING_MPI
	int SizeUsed;
#endif
	static DIR_BOUND_COND_PAR pDirBoundCondPar;
	static INIT_COND_PAR pInitCondPar;
protected:
	SUserInterfacePar User;
	static CTimePartitioning * pTimePar;
	static CLevelPar * pTheHighestLevel;
	void UserBuild();
	void RestrictedUserBuild();
	void UserSetUpRHS(const EL_INFO*);
	void ExtendRestrictedMatrix(int, int);
	CLevelPar* GetTheHighestLevel();
	void DirichletBound();
	static REAL UserDirichletBound(const REAL_D);
	static REAL UserInitialCondition(const REAL_D);
	void UserBuildStatic();
	void Solve();
public:
	void InitialCondition();
	static REAL* GetOldSolution(const REAL_D);
	CLevelPar(int PrevLevel, CUserDataPar *, CLevelPar *);
	void PostprocessMicroSolution();
	void SaveSolutionToVTK(const char*);
	virtual ~CLevelPar();
	void Solve(bool WantAssemble);
	REAL CalcL2Error();
};

#endif // #ifdef CLEVEL_PAR_H
