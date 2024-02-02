#ifndef FEHMM_H
#define FEHMM_H

#ifndef _ALBERTA_H_
#include "alberta.h"
#endif

#include "CLevel.h"
#include "SaveToVTK.h"

/// feHMM is main class encapsulating all the functionality
// feHMM communicates with the user

class feHMM {
protected:
	CLevel* pLevels; ///< the chain list of all levels
public:
	feHMM();
	feHMM(CUserData *); ///<
	~feHMM(); ///<
	void Solve(bool WantAssemble = false);
	void SaveSolutionToVTK(const char*);
	REAL L2Error();
};

#endif // #ifndef FEHMM_H
