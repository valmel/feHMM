#ifndef CLEVEL_H
#define CLEVEL_H

#include "alberta.h"

#include "CUserData.h"
#include "ErrorReturnCodes.h"
#include "SaveToVTK.h"
#include "CAlbertaInterface.h"

#ifdef USING_MPI
#include "mpi.h"
#endif

/// CLevel - class describing one level (macro, mezo, micro...)

//typedef REAL* (*TO_INTERPOLATE)(const FE_SPACE*, DOF_REAL_VEC*, DOF_REAL_VEC*);
//typedef void (*INTERPOLATE)(const FE_SPACE*, DOF_REAL_VEC*, DOF_REAL_VEC*);

typedef REAL (*DIRICHLET_FUNC)(const REAL_D);

struct SUserInterface {
	BUILD_STATIC_MATRIX_TO_PASS pBuildStatic;
	BUILD_DYNAMIC_MATRIX_TO_PASS pBuildDynamic;
	BUILD_RHS_TO_PASS pSetUpRHS;
	MICRO_SOL_POSTPROCESS pMicroSolPostprocess;
    EXACT_SOL pExact;
};

struct SAuxiliaryData {
	int *PrevNumOfBasFct;
	const EL_INFO* PrevElInfo;
	int IndexOfBasFct;
	int IndexOfFEspace;
};

#ifdef USING_MPI
struct SMPIInfo {
	int Rank;
	int Size;
	MPI_Status Status;
};

enum EMPIFlags
{
	EMF_OK = 0, // default status flag
	EMF_CALCULATE, // data send are for calculation
	EMF_STOP // terminate child processes
};

#endif

class CLevel {
	static EXACT_SOL pExact;
	static DIR_BOUND_COND pDirBoundCond;
	static REAL UserLocalExact(const REAL_D);
    SaveToVTK * pSaveVTK;

protected:
	virtual REAL LocalExact(const REAL_D);
	static int ActComponentIdx;

	int IsAssembled;
	const EL_INFO* pElInfo;
	int NumOfLevels; ///< pocet urovni
	int ActualLevel; ///< cislo aktualnej urovne. 0 najnizsia uroven (mikro)

	static CLevel * pPrevLevelStatic;
	CLevel * pPrevLevel;
	CLevel * pNextLevel; ///< smernik na nasledujucu uroven. NULL, ak je posledna

	REAL ScaleFactor;
	REAL L2Error;
	DOF_REAL_VEC ** pSolution;

	SolverInfo Solver;
	CAlbertaInterface* pAlberta;
	SUserInterface User;
	static SAuxiliaryData Aux;

	static int SpaceIndex;

	const int** pMatrixMask; ///< at each level says which parts of ProblemMatrix should ne aasembled - DimOfProb x DimOfProb

#ifdef USING_MPI
	double * pCommunicationMatrix; // = NULL;
	double * pCommunicationVector;
	int MatrixSize;
	SMPIInfo MPIInfo;

	void AssembleMaster();
	void AssembleSlave();

	void SetUpRHSMaster();
	void SetUpRHSSlave();
#endif
	/*members*/
	virtual void Initialize(CUserData *);
	void InterpolateFromPreviousLevel(const EL_INFO*, int, int);
	static REAL Dirichlet(const REAL_D);

	virtual REAL Integrate(int, int, int, int);
	void SolveInt();
	virtual void UserBuild();
	virtual void UserBuildStatic();
	virtual void RestrictedUserBuild();
	virtual void UserSetUpRHS(const EL_INFO *);
	virtual void ExtendRestrictedMatrix(int, int);
	void SetUpRHSForOneElement(const EL_INFO *);
	int GetNumOfLeafEl();
	void ComputeForElInfo(const EL_INFO *);
	void SetUpScaleFactor(REAL);
	static REAL UserDirichletBound(const REAL_D);
	virtual void PostprocessMicroSolution();
	virtual void Solve();
	CLevel() : pSaveVTK(NULL) {
	}
	;
public:
	CLevel* GetNextLevel() {
		return pNextLevel;
	}
	CLevel(int PrevLevel, CUserData *, CLevel *);
	virtual ~CLevel();
	virtual int Assemble(); ///< system matrix assembling
	virtual void Solve(bool WantAssemble);
	virtual void SaveSolutionToVTK(const char*);
	int SetUpRHS();
	virtual void DirichletBound();
	DOF_REAL_VEC** GetSolution();
	virtual REAL CalcL2Error();
	REAL GetL2Error(){
		return L2Error;
	};
	void PrintSolution(int);
};

#endif // #ifdef CLEVEL_H
