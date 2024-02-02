#include "CLevelPar.h"

CLevelPar* CLevelPar::pTheHighestLevel= NULL;

DIR_BOUND_COND_PAR CLevelPar::pDirBoundCondPar= NULL;
INIT_COND_PAR CLevelPar::pInitCondPar= NULL;

EXACT_SOL_PAR CLevelPar::pExact= NULL;
REAL CLevelPar::Time = 0;

REAL** CLevelPar::LocalDofVec;
REAL* CLevelPar::OldSolution;
int CLevelPar::DimOfProb = 0;

CTimePartitioning CLevelPar::TimePar = { 0, 0, 0, 0 };

using namespace std;

CLevelPar::CLevelPar(int PrevLevel/**previous level*/,
		CUserDataPar * Data/**user data*/, CLevelPar *PrevLevelPtr) :
	CLevel() {

	ActualLevel = PrevLevel-1;
	pPrevLevel = PrevLevelPtr;

	NumOfLevels = Data->NumOfLevels;
	ScaleFactor = Data->pScaleFactor[NumOfLevels-ActualLevel-1];
	DimOfProb = Data->DimOfProb;

	pAlberta = new CAlbertaInterfacePar();

	IsAssembled = false;

	User.pBuildStatic = NULL;
	User.pBuildDynamic = NULL;
	User.pSetUpRHS = NULL;
	User.pExact = NULL;

	pSaveVTKpar = NULL;
	IsBiLinFormDependentOnPrevStep=Data->IsBiLinFormDependentOnPrevStep;

	if (ActualLevel==0) {
		User.pBuildStatic = Data->pBuildStatic;
		User.pBuildDynamic = Data->pBuildDynamic;
		User.pMicroSolPostprocess = Data->pMicroSolPostprocess;
		pTheHighestLevel = GetTheHighestLevel();
	}
	if (ActualLevel==NumOfLevels-1) {
		TimePar.StartTime = Data->StartTime;
		TimePar.EndTime = Data->EndTime;
		TimePar.Tau = Data->TimeStep;

		pExact = Data->pExactSolution;
		User.pExact = Data->pExactSolution;
		pSaveVTKpar = new SaveToVTKpar(Data->pExactSolution, Data->DimOfWorld, Data->DimOfProb);

		User.pSetUpRHS = Data->pSetUpRHS;
		pDirBoundCondPar = Data->pDirBoundCond;
		pInitCondPar = Data->pInitialCondition;
	}

	Initialize(Data);

	if (ActualLevel)
		pNextLevel = new CLevelPar(ActualLevel, Data, this);

	if (!ActualLevel)
		pNextLevel = NULL;

	LocalDofVec = (REAL**)malloc(pAlberta->DimOfProb*sizeof(REAL*));
	for (int i = 0; i < pAlberta->DimOfProb; i++)
		LocalDofVec[i] = (REAL*)malloc(pAlberta->NumOfBasFct[i]*sizeof(REAL));

	OldSolution = (REAL*)malloc(pAlberta->DimOfProb*sizeof(REAL));

	pAlberta->FillRestrictedVertexIndex(); //robi to aj v prvom leveli a co ? ;)


	if (ActualLevel==NumOfLevels-1) {
		SetUpScaleFactor(1.);
		/* Check whether to make static build */
		CLevelPar *pLevel;
		pLevel=this;
		while (pLevel->pNextLevel)
			pLevel=(CLevelPar*)pLevel->pNextLevel;
		pLevel->UserBuildStatic();
		/*
		 pAlberta->CreateStaticMatrices();
		 printf("Static build\n");
		 UserBuildStatic();
		 */
	}

#ifdef USING_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIInfo.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIInfo.Size);

	MatrixSize = 0;
	SizeUsed = 0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		MatrixSize += pAlberta->NumOfBasFct[i];
		SizeUsed += pAlberta->FEspace[i]->admin->size_used;
	}
	pCommunicationMatrix = NULL;
	pCommunicationMatrix
			= (double*)malloc(MatrixSize*MatrixSize*sizeof(double));
	pCommunicationVector = NULL;
	pCommunicationVector = (double*)malloc(MatrixSize*sizeof(double));
#endif
}

CLevelPar::~CLevelPar() {
#ifdef USING_MPI
	if (pCommunicationMatrix) {
		free(pCommunicationMatrix);
		pCommunicationMatrix = NULL;
	}
	if (pCommunicationVector) {
		free(pCommunicationVector);
		pCommunicationVector = NULL;
	}
#endif

	if (pSaveVTKpar) {
		delete pSaveVTKpar;
		pSaveVTKpar=NULL;
	}

	if (LocalDofVec) {
		for (int i = 0; i < pAlberta->DimOfProb; i++)
			free(LocalDofVec[i]);
		free(LocalDofVec);
		LocalDofVec = NULL;
	}

	if (OldSolution) {
		free(OldSolution);
		OldSolution = NULL;
	}
}

void CLevelPar::UserBuild() {
	DOF_REAL_VEC** vec = ((CAlbertaInterfacePar*)(pTheHighestLevel->pAlberta))->UhOld;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j=0; j < vec[i]->fe_space->bas_fcts->n_bas_fcts; j++) {
			LocalDofVec[i][j] = vec[i]->vec[(pTheHighestLevel->pAlberta->Dofs[i])[j]];
		}
	}
	if (User.pBuildDynamic) {
		User.pBuildDynamic(pAlberta->FEspace, pAlberta->Matrix,
				&CLevelPar::GetOldSolution, TimePar.Time, TimePar.Tau);
		pAlberta->SumMatrices(pMatrixMask);
	}
	if (!User.pBuildDynamic && User.pBuildStatic) { // second should be obvious
		pAlberta->AssembleProblemMatrix(pMatrixMask,
				(const DOF_MATRIX ***)pAlberta->StaticMatrix);
	} else
		pAlberta->AssembleProblemMatrix(pMatrixMask,
				(const DOF_MATRIX ***)pAlberta->Matrix);
}

void CLevelPar::RestrictedUserBuild() {
	if (User.pBuildDynamic) {
		User.pBuildDynamic(pAlberta->RestrictedFEspace,
				pAlberta->RestrictedMatrix, &CLevelPar::GetOldSolution,
				TimePar.Time, TimePar.Tau);
		pAlberta->SumRestrictedMatrices(pMatrixMask);
	}
}

void CLevelPar::ExtendRestrictedMatrix(int i, int j) {
	if (ActualLevel == 0 && !User.pBuildDynamic)
		pAlberta->ExtendRestrictedMatrix(i, j,
				(const DOF_MATRIX ***)pAlberta->StaticRestrictedMatrix);
	else
		pAlberta->ExtendRestrictedMatrix(i, j,
				(const DOF_MATRIX ***)pAlberta->RestrictedMatrix);
}

void CLevelPar::Solve() {
	int idx, i, j;

	idx=0;
	for (i = 0; i < pAlberta->DimOfProb; i++) {
		for (j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->ProblemRHS->vec[idx]=pAlberta->Fh[i]->vec[j];
			idx++;
		}
	}

	oem_solve_s(pAlberta->ProblemMatrix, pAlberta->ProblemRHS,
			pAlberta->ProblemSol, Solver.solver, Solver.tol, Solver.icon,
			Solver.restart, Solver.miter, Solver.info);
	
	idx = 0;
	for (i = 0; i < pAlberta->DimOfProb; i++) {
		for (j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->Uh[i]->vec[j] = pAlberta->ProblemSol->vec[idx];
			((CAlbertaInterfacePar*)pAlberta)->UhOld[i]->vec[j]
					= pAlberta->ProblemSol->vec[idx];
			idx++;
		}
	}

#ifdef USING_MPI	
	MPI_Bcast(pAlberta->ProblemSol->vec, SizeUsed, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (MPIInfo.Rank != 0) {
		idx=0;
		for (i = 0; i < pAlberta->DimOfProb; i++) {
			for (j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
				pAlberta->Uh[i]->vec[j] = pAlberta->ProblemSol->vec[idx];
				((CAlbertaInterfacePar*)pAlberta)->UhOld[i]->vec[j]
						= pAlberta->ProblemSol->vec[idx];
				idx++;
			}
		}
	}
#endif
}

CLevelPar* CLevelPar::GetTheHighestLevel() {
	if (pPrevLevel)
		return ((CLevelPar*)pPrevLevel)->GetTheHighestLevel();
	else
		return this;
}

REAL* CLevelPar::GetOldSolution(const REAL_D x) {
	const EL_INFO * el_info = pTheHighestLevel->pElInfo;
	DOF_REAL_VEC** vec = ((CAlbertaInterfacePar*)(pTheHighestLevel->pAlberta))->UhOld;
	REAL bary[N_LAMBDA];
	int is = 0;

	is = world_to_coord(el_info, x, bary);
	if (is != -1)
		exit(0);

	for (int i = 0; i < DimOfProb; i++)
		OldSolution[i] = eval_uh(bary, LocalDofVec[i],
				vec[i]->fe_space->bas_fcts);

	return OldSolution;
}

void CLevelPar::UserSetUpRHS(const EL_INFO *pElInfo) {
	DOF_REAL_VEC** vec = ((CAlbertaInterfacePar*)(pTheHighestLevel->pAlberta))->UhOld;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j=0; j < vec[i]->fe_space->bas_fcts->n_bas_fcts; j++) {
			LocalDofVec[i][j]=vec[i]->vec[(pTheHighestLevel->pAlberta->Dofs[i])[j]];
		}
	}
	User.pSetUpRHS(pElInfo, pAlberta->FEspace, pAlberta->LocalFh,
			&CLevelPar::GetOldSolution, TimePar.Time, TimePar.Tau);
}

REAL CLevelPar::UserDirichletBound(const REAL_D x) {
	return pDirBoundCondPar(x, TimePar.Time, SpaceIndex);
}

void CLevelPar::DirichletBound() {
	if (pDirBoundCondPar)
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			SpaceIndex = i;
			dirichlet_bound(&CLevelPar::UserDirichletBound, pAlberta->Fh[i],
					pAlberta->Uh[i], nil);
		}
}

REAL CLevelPar::UserInitialCondition(const REAL_D x) {
	return pInitCondPar(x, TimePar.Time, SpaceIndex);
}

void CLevelPar::InitialCondition() {
	if (pInitCondPar)
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			SpaceIndex = i;
			interpol(&CLevelPar::UserInitialCondition, pAlberta->Uh[i]);
			interpol(&CLevelPar::UserInitialCondition, ((CAlbertaInterfacePar*)pAlberta)->UhOld[i]);
		}

	int idx = 0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->ProblemSol->vec[idx] = pAlberta->Uh[i]->vec[j];
			idx++;
		}
	}
}

void CLevelPar::PostprocessMicroSolution() {
	if (!User.pMicroSolPostprocess)
		return;

	int idx=0;

	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->Uh[i]->vec[j]
					= pAlberta->ProblemUh[Aux.IndexOfFEspace][Aux.IndexOfBasFct]->vec[idx];
			idx++;
		}
	}
	User.pMicroSolPostprocess(pAlberta->FEspace, pAlberta->Uh, TimePar.Time,
			TimePar.Tau);
}

void CLevelPar::UserBuildStatic() {
	REAL_D x;

	if (!User.pBuildStatic)
		return;
	if (!pAlberta->StaticMatrix)
		pAlberta->CreateStaticMatrices();
	for (int i = 0; i < pAlberta->DimOfWorld; i++)
		x[i] = (pAlberta->x[i] + pAlberta->y[i])/2.+0.01;
	pAlberta->TransformMeshes(x);
	User.pBuildStatic(pAlberta->FEspace, pAlberta->StaticMatrix,
			&CLevelPar::GetOldSolution, TimePar.Time, TimePar.Tau);
	User.pBuildStatic(pAlberta->RestrictedFEspace,
			pAlberta->StaticRestrictedMatrix, &CLevelPar::GetOldSolution,
			TimePar.Time, TimePar.Tau);
}

void CLevelPar::Solve(bool WantAssemble) {
	FUNCNAME("CLevelPar::Solve(bool WantAssemble)");

	L2Error = 0;

#ifdef USING_MPI
	if (MPIInfo.Rank == 0)
#endif
		INFO(9,9,"Solve\n");

	TimePar.Time = TimePar.StartTime;
#ifdef USING_MPI
	if (MPIInfo.Rank == 0)
#endif
		InitialCondition(); // initial condition must be first...
	
#ifdef USING_MPI
	int idx;
	MPI_Bcast(pAlberta->ProblemSol->vec, SizeUsed, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (MPIInfo.Rank != 0) {
		idx=0;
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
				pAlberta->Uh[i]->vec[j] = pAlberta->ProblemSol->vec[idx];
				((CAlbertaInterfacePar*)pAlberta)->UhOld[i]->vec[j]
						=pAlberta->ProblemSol->vec[idx];
				idx++;
			}
		}
	}
#endif
	
	if (!IsAssembled || WantAssemble) {
#ifdef USING_MPI
		if (MPIInfo.Rank == 0)
#endif
			if (!IsAssembled)
				INFO(9,9,"First assembling of system matrix\n");
		Assemble();
		IsAssembled = true;
	}

	TimePar++;

#ifdef USING_MPI
	if (MPIInfo.Rank == 0)
#endif
		printf("Passed initial condition\n");

	for (; TimePar.Time <= TimePar.EndTime; TimePar++) {
		if (IsBiLinFormDependentOnPrevStep)
			if (TimePar.Time != TimePar.StartTime + TimePar.Tau)
				Assemble();
		SetUpRHS();
		DirichletBound();
		Solve();
		L2Error += TimePar.Tau * pow(CalcL2Error(),2.);
	}
	L2Error = pow(L2Error, 0.5);
	TimePar--;
}

void CLevelPar::SaveSolutionToVTK(const char *pFileName) {
	pSolution = GetSolution();
	pSaveVTKpar->SaveVTKfile((const DOF_REAL_VEC **)pSolution, pFileName, TimePar.Time);
}

REAL CLevelPar::UserLocalExact(const REAL_D x) {
	return pExact(x, Time, ActComponentIdx);
}

REAL CLevelPar::LocalExact(const REAL_D x) {
	return UserLocalExact(x);
}

REAL CLevelPar::CalcL2Error() {
	FUNCNAME("CLevelPar::CalcL2Error()");
	
	REAL RetVal, Err;
	pSolution = GetSolution();
	RetVal = 0;
	Time = TimePar.Time;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		ActComponentIdx = i;
		Err = L2_err(&CLevelPar::UserLocalExact, pSolution[i], nil, true, nil,
				nil);
		RetVal += Err*Err;
	}
	
#ifdef USING_MPI
		if (MPIInfo.Rank == 0)
#endif
		INFO(9,9,"L2error local = %f\n",pow(RetVal, 0.5));
	return pow(RetVal, 0.5);
}
