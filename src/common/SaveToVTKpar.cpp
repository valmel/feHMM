#include "SaveToVTKpar.h"

using namespace std;

REAL SaveToVTKpar::GlobalTime = 0;

REAL SaveToVTKpar::LocalExact(const REAL_D x) {
	return (pExact) ? pExact(x, GlobalTime, ActComponentIdx) : 0;
}

SaveToVTKpar::SaveToVTKpar(EXACT_SOL_PAR exact, int world, int prob) : SaveToVTK((EXACT_SOL)1, world, prob) {
	if (exact)
		pExact=exact;
}

void SaveToVTKpar::SaveVTKfile(const DOF_REAL_VEC **pSol, const char *pName, REAL Time) {
	GlobalTime = Time;
	SaveToVTK::SaveVTKfile(pSol, pName);
}

