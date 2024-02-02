#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "feHMM.h"

using namespace std;

feHMM::feHMM() {
}

/**
   Constructor checks user data and initializes CLevel array
*/
feHMM::feHMM(CUserData * Data) {
	Data->CheckUserData(); // terminates program if failed
	pLevels = new CLevel(Data->NumOfLevels, Data, NULL);
}

/**
   Destructor goes throught all levels and destructs them
*/
feHMM::~feHMM() {
	CLevel* temp = pLevels;
	while (temp) {
		pLevels = pLevels->GetNextLevel();
		delete temp;
		temp = pLevels;
	}
	pLevels = NULL;
}

void feHMM::Solve(bool WantAssemble) {
	if (pLevels)
		pLevels->Solve(WantAssemble);
	else { // this realy should not happened
		ERROR("pLevels not initialized. Possible problem in the library.\n");
		exit(ERC_LEVELS_NOT_INITIALIZED);
	}
}

void feHMM::SaveSolutionToVTK(const char * pName) {
    if (pLevels)
    	pLevels->SaveSolutionToVTK(pName);
	else { // this realy should not happened
		ERROR("pLevels not initialized. Possible problem in the library.\n");
		exit(ERC_LEVELS_NOT_INITIALIZED);
	}
}

REAL feHMM::L2Error(){
	REAL RetVal;
	if (pLevels)
	    	RetVal = pLevels->GetL2Error();
		else { // this realy should not happened
			ERROR("pLevels not initialized. Possible problem in the library.\n");
			exit(ERC_LEVELS_NOT_INITIALIZED);
	}
	return RetVal;
}
