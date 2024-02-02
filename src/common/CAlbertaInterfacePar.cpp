#include "CAlbertaInterfacePar.h"

using namespace std;

void CAlbertaInterfacePar::Initialize(
		CUserData * Data/** Data provided by user */, int ActLevel,
		int Rank) {
	CAlbertaInterface::Initialize(Data, ActLevel, Rank);
	UhOld = NULL;
	if (ActualLevel == NumOfLevels-1) {
		UhOld = (DOF_REAL_VEC**)malloc(DimOfProb
				*sizeof(DOF_REAL_VEC*));
		for (int i = 0; i < DimOfProb; i++)
			UhOld[i] = get_dof_real_vec("Uh", FEspace[i]);
	}
}

void CAlbertaInterfacePar::Free() {
	if (UhOld) {
		for (int i = 0; i < DimOfProb; i++)
			free_dof_real_vec(UhOld[i]);
		free(UhOld);
		UhOld = NULL;
	}
	CAlbertaInterface::Free();
}

