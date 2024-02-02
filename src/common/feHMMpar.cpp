#include "feHMMpar.h" 

using namespace std;

/**
 Konstruktor skontroluje uzivatelske data a vytvori zretazeny zoznam jednotlivych
 urovni (najvyssiu uroven, zvysne sa vytvoria samostatne).
 */

feHMMpar::feHMMpar(CUserDataPar* Data) :
	feHMM() {
	Data->CheckUserData();  // terminates program if failed
	pLevels = new CLevelPar(Data->NumOfLevels, Data, NULL);
}

/**
 Destruktor prechadzanim do hlbky prejde cez vsetky urovne a zavola na nich delete.
 */

feHMMpar::~feHMMpar() {
	CLevel* temp = pLevels;
	while (temp) {
		pLevels = pLevels->GetNextLevel();
		delete temp;
		temp = pLevels;
	}
	pLevels = NULL;
}
