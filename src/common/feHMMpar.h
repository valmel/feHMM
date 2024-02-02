#ifndef FEHMM_PAR_H
#define FEHMM_PAR_H

#include "feHMM.h"
#include "CLevelPar.h"
#include "CTimePartitioning.h"
#include "CUserDataPar.h"

class feHMMpar : public feHMM {
public:
	feHMMpar(CUserDataPar *);
	~feHMMpar();
};

#endif // #ifdef FEHMM_PAR_H
