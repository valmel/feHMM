#ifndef CTIME_PARTITIONING_H
#define CTIME_PARTITIONING_H
#include "alberta.h"

class CTimePartitioning{
public:
	REAL StartTime; ///>the beginning of time interval
	REAL EndTime; ///>the end of time interval
	REAL Tau; ///> the time step
	REAL Time; ///> actual time
	const CTimePartitioning& operator++(int){Time+=Tau; return *this;};
	const CTimePartitioning& operator--(int){Time-=Tau; return *this;};
};

#endif //#ifndef CTIME_PARTITIONING_H

