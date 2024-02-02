#ifndef CALBERTA_INTERFACE_PAR_H
#define CALBERTA_INTERFACE_PAR_H

#include "CAlbertaInterface.h"
#include "CUserDataPar.h"

class CAlbertaInterfacePar: public CAlbertaInterface {
public:
	DOF_REAL_VEC** UhOld;
 	void Initialize(CUserData *, int, int);
	void Free();
};

#endif // #ifndef CALBERTA_INTERFACE_PAR_H

