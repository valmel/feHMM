#ifndef SAVE_TO_VTK_PAR_H
#define SAVE_TO_VTK_PAR_H

#include "SaveToVTK.h"
#include "CUserDataPar.h"

class SaveToVTKpar : public SaveToVTK{
	EXACT_SOL_PAR pExact;
	REAL LocalExact(const REAL_D);
	static REAL GlobalTime;
public:
	void SaveVTKfile(REAL);
	SaveToVTKpar(EXACT_SOL_PAR, int, int);
	virtual void SaveVTKfile(const DOF_REAL_VEC **, const char *, REAL);
};

#endif // #ifndef SAVE_TO_VTK_PAR_H
