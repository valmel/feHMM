#include "CUserData.h"
#include <stdio.h>
#include <iostream>
#include <cstdarg>

using namespace std;

CUserData::~CUserData() {
	if (pNumOfRefin)
		delete []pNumOfRefin;
	if (pQuadDeg)
		delete []pQuadDeg;
	if (pScaleFactor)
		delete []pScaleFactor;
	if (pFEdeg)
		delete []pFEdeg;
	if (pMacroFile)
		delete []pMacroFile;
}

void CUserData::ExtractParamFromLine(const char * pLine, int nParam,
		int Type, void * pDst) {
	double *pDouble;
	float FloatVal;
	int *pInt;
	char TmpStr[100], *pStr;
	int Idx, Jdx;
	pDouble=NULL;
	pInt=NULL;
	switch (Type) {
	case 1:
		pInt=(int*)pDst;
		break;
	case 2:
		pDouble=(double*)pDst;
		break;
	}
	pStr=(char*)pLine;
	for (Idx=0; Idx<nParam; Idx++) {
		Jdx=0;
		while (!isspace(*pStr) && (*pStr!='\0') && (*pStr!=',')) {
			TmpStr[Jdx++]=*pStr;
			pStr++;
		}
		pStr++; // skip ','
		TmpStr[Jdx++]='\0';
		switch (Type) {
		case 1:
			sscanf(TmpStr, "%d", &pInt[Idx]);
			break;
		case 2:
			sscanf(TmpStr, "%f", &FloatVal);
			pDouble[Idx]=(double)FloatVal;
			break;
		}
	}
}

void CUserData::LoadUserDataFromFile(const char * pFileName) {
	FUNCNAME("CUserData::LoadUserDataFromFile");
	init_parameters(0, pFileName);
	char TmpStr[1000];

	GET_PARAMETER(1, "number of levels", "%d", &NumOfLevels);
	GET_PARAMETER(1, "dim of world", "%d", &DimOfWorld);
	GET_PARAMETER(1, "dim of problem", "%d", &DimOfProb);
	GET_PARAMETER(1, "problem number", "%d", &ProblemNum);
	GET_PARAMETER(1, "boundary condition", "%d", &BoundCond);

	pQuadDeg = new int[NumOfLevels];
	pScaleFactor = new double[NumOfLevels];
	pNumOfRefin = new int[NumOfLevels];
	pFEdeg = new int[DimOfProb];

	GET_PARAMETER(1, "quadrature", "%s", TmpStr);
	ExtractParamFromLine(TmpStr, NumOfLevels, 1, pQuadDeg);
	GET_PARAMETER(1, "scaling factors", "%s", TmpStr);
	ExtractParamFromLine(TmpStr, NumOfLevels, 2, pScaleFactor);
	GET_PARAMETER(1, "level refinement", "%s", TmpStr);
	ExtractParamFromLine(TmpStr, NumOfLevels, 1, pNumOfRefin);
	GET_PARAMETER(1, "finite element degree", "%s", TmpStr);
	ExtractParamFromLine(TmpStr, DimOfProb, 1, pFEdeg);
	GET_PARAMETER(1, "macro file", "%s", TmpStr);

	pMacroFile = new char[strlen(TmpStr)+1];
	sprintf(pMacroFile, "%s", TmpStr);
	// solver
	GET_PARAMETER(1, "solver", "%d", &Solver.solver);
	GET_PARAMETER(1, "solver tolerance", "%f", &Solver.tol);
	GET_PARAMETER(1, "solver precon", "%d", &Solver.icon);
	GET_PARAMETER(1, "solver max iteration", "%d", &Solver.miter);
	GET_PARAMETER(1, "solver info", "%d", &Solver.info);
	if (Solver.solver == GMRes)
		GET_PARAMETER(1, "solver restart", "%d", &Solver.restart);

#ifdef USING_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#endif
}

int CUserData::CheckUserData() {
	FUNCNAME("CUserData::CheckUserData");
	
	int RetVal;
	
	RetVal = CheckUserDataWithoutFuncPointers();
	if (RetVal != ERC_OK)
		exit(RetVal);
	
	if (!pBuildStatic && !pBuildDynamic) {
		ERROR("Building function is not set. At least one of CUserData.pBuildStatic, CUserData.pBuildDynamic has to be set.\n");
		exit(ERC_BUILD_FUNC_MISSING);
	}
	if (!pSetUpRHS) {
		ERROR("CUserData.pSetUpRHS is not set.\n");
		exit(ERC_SETUPRHS_FUNC_MISSING);
	}
	if (!pExactSolution)
		WARNING("CUserData.pExactSolution is not set.\n");
	if (!pDirBoundCond)
		WARNING("CUserData.pDirBoundCond is not set.\n");
	return ERC_OK;
}

int CUserData::CheckUserDataWithoutFuncPointers() {
	FUNCNAME("CUserData::CheckUserDataWithoutFuncPointers");
	if (NumOfLevels<0) {
		ERROR("CUserData.NumOfLevels is not set.\n");
		exit(ERC_NUM_OF_LEVELS_MISSING);
	}
	
	if (DimOfWorld < 1 || DimOfWorld > 3 ) {
		ERROR("CUserData.DimOfWorld is not correctly set (is not one of {1,2,3}).\n");
		exit(ERC_DIM_OF_WORLD_MISSING);
	}
	if (DimOfProb < 1) {
		ERROR("CUserData.DimOfProb is not set.\n");
		exit(ERC_DIM_OF_PROB_MISSING);
	}
	if (!pNumOfRefin) {
		ERROR("CUserData.pNumOfRefin is not set.\n");
		exit(ERC_REFINEMENT_MISSING);
	}
	if (!pQuadDeg) {
		ERROR("CUserData.pQuadDeg is not set.\n");
		exit(ERC_QUADRATURE_MISSING);
	}
	if (!pScaleFactor) {
		ERROR("CUserData.pScaleFactor is not set.\n");
		exit(ERC_SCALE_FACTOR_MISSING);
	}
	if (!pFEdeg) {
		ERROR("CUserData.pFEdeg is not set.\n");
		exit(ERC_FE_DEG_MISSING);
	}
	if (!pMacroFile) {
		ERROR("CUserData.MacroFile is not set.\n");
		exit(ERC_ERROR_MACRO_FILE_MISSING);
	}
	if (!Solver.solver) {
		ERROR("CUserData.Solver.solver is not set.\n");
		exit(ERC_SOLVER_NOT_SET);
	}
	if (!pMatrixMask) {
		ERROR("MatrixMask is not set.\n");
		exit(ERC_MATRIX_MASK_NOT_SET);
	}
	if (BoundCond < 1) {
		WARNING("Boundary condition for microproblems is not set. The periodic one is used.\n");
		BoundCond=1;
	}
	return ERC_OK;
}

void CUserData::FillSolverData(SolverInfo & Dst) {
	Dst.solver = Solver.solver;
	Dst.tol = Solver.tol;
	Dst.icon = Solver.icon;
	Dst.miter = Solver.miter;
	Dst.info = Solver.info;
	if (Dst.solver == GMRes)
		Dst.restart = Solver.restart;
}
