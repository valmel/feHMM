#include "CUserDataPar.h"

#define RANDOM_NUMBER  -123454321

using namespace std;

void CUserDataPar::LoadUserDataFromFile(const char * pFileName) {

	CUserData::LoadUserDataFromFile(pFileName);

	REAL Val = 0;

	StartTime = EndTime = TimeStep = NumTimeStep = RANDOM_NUMBER;
	if (GET_PARAMETER(1, "start time", "%f", &Val))
		StartTime = Val;
	if (GET_PARAMETER(1, "end time", "%f", &Val))
		EndTime = Val;
	if (GET_PARAMETER(1, "time step", "%f", &Val))
		TimeStep = Val;
	GET_PARAMETER(1, "number of time steps", "%d", &NumTimeStep);
}

int CUserDataPar::CheckUserData() {
	FUNCNAME("CUserDataPar::CheckUserData");

	int RetVal = CheckUserDataWithoutFuncPointers();

	if (RetVal != ERC_OK)
		exit(RetVal);

	if (!pBuildStatic && !pBuildDynamic) {
		ERROR("At least of CUserDataPar.pBuildStatic, CUserDataPar.pBuildDynamic has to be set.\n");
		exit(ERC_BUILD_FUNC_MISSING);
	}
	if (!pSetUpRHS) {
		ERROR("CUserDataPar.pSetUpRHS is not set.\n");
		exit(ERC_SETUPRHS_FUNC_MISSING);
	}
	if (!pExactSolution)
		WARNING("CUserDataPar.pExactSolution is not set.\n");
	if (!pDirBoundCond)
		WARNING("CUserDataPar.pDirBoundCond is not set.\n");
	if (!pInitialCondition)
		WARNING("CUserDataPar.pInitialCondition is not set. Starting from zeros.\n");
	if (IsBiLinFormDependentOnPrevStep < 0) {
		ERROR("CUserDataPar.IsBiLinFormDependentOnPrevStep is not set.\n");
		exit(ERC_IS_BI_LIN_FORM_DEPENDENT_ON_PREV_STEP_MISSING);
	}
	if (StartTime == RANDOM_NUMBER || NumTimeStep == RANDOM_NUMBER){
		if (StartTime == RANDOM_NUMBER)
			ERROR("Start time is not set.\n");
		else
			ERROR("Number of time steps is not set.\n");
		exit(ERC_TIME_INFO_NOT_COMPLETE);
	}
	if (EndTime == RANDOM_NUMBER && TimeStep == RANDOM_NUMBER){
		ERROR("Both EndTime and TimeStep are not set. At least one of them has to be set.\n");
		exit(ERC_TIME_INFO_NOT_COMPLETE);
	}
	if (NumTimeStep == 0){
		ERROR("Number of time steps is zero\n");
		exit(ERC_ZERO_NUM_OF_TIME_STEP);
	}
	if (TimeStep != RANDOM_NUMBER && TimeStep!=0)
		EndTime = StartTime + TimeStep * NumTimeStep;
	if (EndTime != RANDOM_NUMBER && EndTime != StartTime)
		TimeStep = (EndTime - StartTime)/(1.0 * NumTimeStep);
	if (EndTime == StartTime && TimeStep == 0){
		ERROR("StartTime equals EndTime and TimeStep is zero.\n");
		exit(ERC_TIME_INFO_NOT_COMPLETE);
	}

#ifdef USING_MPI
	if (Rank != 0) // write info only on master
		return ERC_OK;
#endif
	if (EndTime != RANDOM_NUMBER && TimeStep != RANDOM_NUMBER && TimeStep != 0)
		INFO(9,9,"Both TimeStep and EndTime are set. EndTime is overriden by TimeStep*NumTimeStep\n");

	INFO(1,1,"Start time: %f\n",StartTime);
	INFO(1,1,"  End time: %f\n",EndTime);
	INFO(1,1," Time step: %f\n",TimeStep);
	return ERC_OK;
}

