#ifndef SAVE_TO_VTK_H
#define SAVE_TO_VTK_H

#include <stdio.h>
#include <iostream>
#include "alberta.h"
#include "CUserData.h"
#include "ErrorReturnCodes.h"

class SaveToVTK {
	EXACT_SOL pExact;

	int DimOfProb;
	int DimOfWorld;
	DOF_INT_VEC * pDofInd;
	MESH * pMesh;
	MACRO_DATA * pMacroData;

	void FillDofIndices(DOF_INT_VEC *, MESH *);
protected:
	static int ActComponentIdx;
	virtual REAL LocalExact(const REAL_D);
public:
	SaveToVTK(EXACT_SOL, int, int);
	virtual ~SaveToVTK();
	virtual void SaveVTKfile(const DOF_REAL_VEC **, const char *);
};

#endif // #ifndef  SAVE_TO_VTK_H
