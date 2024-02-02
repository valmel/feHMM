#include <stdio.h>
#include <stdlib.h>
#include "CLevel.h"
#include <iostream>

using namespace std;

SAuxiliaryData CLevel::Aux = { NULL, NULL, 0, 0 };

CLevel* CLevel::pPrevLevelStatic= NULL;
EXACT_SOL CLevel::pExact = NULL;
int CLevel::ActComponentIdx = 0;

DIR_BOUND_COND CLevel::pDirBoundCond= NULL;

int CLevel::SpaceIndex = 0;

#ifdef USING_MPI
#include "CLevelMPI.cpp"
#endif // #ifdef USING_MPI
/**
 Level data are set and next level is created if needed. Previous level and data
 structure are on input
 */
CLevel::CLevel(int PrevLevel/**previous level*/,
		CUserData * Data/**user data*/, CLevel *PrevLevelPtr) {

	ActualLevel = PrevLevel-1;
	pPrevLevel = PrevLevelPtr;

	NumOfLevels = Data->NumOfLevels;
	ScaleFactor = Data->pScaleFactor[NumOfLevels-ActualLevel-1];

	pAlberta = new CAlbertaInterface();

	User.pBuildStatic = NULL;
	User.pBuildDynamic = NULL;
	User.pSetUpRHS = NULL;
	User.pExact = NULL;

	pSaveVTK = NULL;

	if (ActualLevel==0) {
		User.pBuildStatic = Data->pBuildStatic;
		User.pBuildDynamic = Data->pBuildDynamic;
		User.pMicroSolPostprocess = Data->pMicroSolPostprocess;
	}
	if (ActualLevel==NumOfLevels-1) {
		pExact = Data->pExactSolution;
		IsAssembled = false;
		User.pSetUpRHS = Data->pSetUpRHS;
		pDirBoundCond = Data->pDirBoundCond;
		User.pExact = Data->pExactSolution;
		pSaveVTK = new SaveToVTK(Data->pExactSolution, Data->DimOfWorld, Data->DimOfProb);
	}

	Initialize(Data);

	if (ActualLevel)
		pNextLevel = new CLevel(ActualLevel,Data,this);

	if (!ActualLevel)
		pNextLevel = NULL;

	if (pNextLevel) {
		pNextLevel->pAlberta->FillRestrictedVertexIndex();
	}

	if (ActualLevel==NumOfLevels-1) {
		SetUpScaleFactor(1.);

		CLevel *pLevel;
		pLevel=this;
		while (pLevel->pNextLevel)
			pLevel=pLevel->pNextLevel;
		pLevel->UserBuildStatic();
	}
	/*
	 if (ActualLevel == 0 && User.pBuildStatic) {
	 pAlberta->CreateStaticMatrices();
	 UserBuildStatic();
	 }
	 */
#ifdef USING_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIInfo.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIInfo.Size);

	MatrixSize = 0;
	for (int i = 0; i < pAlberta->DimOfProb; i++)
		MatrixSize += pAlberta->NumOfBasFct[i];
	pCommunicationMatrix = NULL;
	pCommunicationMatrix
			= (double*)malloc(MatrixSize*MatrixSize*sizeof(double));
	pCommunicationVector = NULL;
	pCommunicationVector = (double*)malloc(MatrixSize*sizeof(double));
#endif
}

/**
 Free all allocated structures
 */
CLevel::~CLevel() {
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

	if (pSaveVTK) {
		delete pSaveVTK;
		pSaveVTK=NULL;
	}

	if (!ActualLevel)
		pNextLevel = NULL;
	if (pAlberta->ProblemUh) {
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			for (int j = 0; j < Aux.PrevNumOfBasFct[i]; j++)
				free_dof_real_vec(pAlberta->ProblemUh[i][j]);
			free(pAlberta->ProblemUh[i]);
		}
		free(pAlberta->ProblemUh);
		pAlberta->ProblemUh = NULL;
	}
	if (ActualLevel) {
		Aux.PrevNumOfBasFct=(int*)realloc(Aux.PrevNumOfBasFct,
				pAlberta->DimOfProb*sizeof(int));
		memcpy(Aux.PrevNumOfBasFct, pAlberta->NumOfBasFct, pAlberta->DimOfProb
				*sizeof(int));
	} else if (Aux.PrevNumOfBasFct) {
		free(Aux.PrevNumOfBasFct);
		Aux.PrevNumOfBasFct=NULL;
	}
	pAlberta->Free(); /* may not be before memcpy above */
	if (pAlberta) {
		delete pAlberta;
		pAlberta = NULL;
	}
}

/**
 Mainly fill pAlberta structure, create mesh and allocate space if needed
 */
void CLevel::Initialize(CUserData * Data/** Data provided by user */) {

#ifdef USING_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIInfo.Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIInfo.Size);

	pAlberta->Initialize(Data, ActualLevel, MPIInfo.Rank);
#else
	pAlberta->Initialize(Data, ActualLevel, 0);
#endif

	pAlberta->ProblemUh = NULL;
	L2Error = 0;

	if (pPrevLevel) {
		pAlberta->ProblemUh = (DOF_REAL_VEC***)malloc(pAlberta->DimOfProb
				*sizeof(DOF_REAL_VEC**));
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			pAlberta->ProblemUh[i]
					= (DOF_REAL_VEC**)malloc(pPrevLevel->pAlberta->NumOfBasFct[i]
							*sizeof(DOF_REAL_VEC*));
			for (int j = 0; j < pPrevLevel->pAlberta->NumOfBasFct[i]; j++) {
				pAlberta->ProblemUh[i][j] = get_dof_real_vec("ProblemUh",
						pAlberta->ProblemFEspace);
			}
		}
	}

	pMatrixMask = (const int**)Data->pMatrixMask;
	Data->FillSolverData(Solver);
}

/**
 On every level (except lowest) lower levels are assembled and solve si called
 afterwards. On lowest level user build is used.
 */
int CLevel::Assemble() {

#ifdef USING_MPI
	if (ActualLevel==NumOfLevels-1) {
		if (MPIInfo.Rank==0)
			AssembleMaster();
		else
			AssembleSlave();
		return 0;
	}
#endif

	TRAVERSE_STACK *stack;
	static REAL Vol;
	static bool RestrictedElement;
	int NumOfLeafEl = 0;
	int iter = 0;
	double pom = 0;

	NumOfLeafEl = GetNumOfLeafEl();
	double OneTenth = NumOfLeafEl/10.0;
	//printf("NUMBER OF LEAF ELEMENTS %d\n", NumOfLeafEl);
	for (int i = 0; i < pAlberta->DimOfProb; i++)
		for (int j = 0; j < pAlberta->DimOfProb; j++) {
			clear_dof_matrix(pAlberta->Matrix[i][j]); //clear dof matrix in the beginning
			clear_dof_matrix(pAlberta->RestrictedMatrix[i][j]); //clear dof matrix in the beginning
		}

	stack = get_traverse_stack();
	pom = OneTenth;
	for (pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL
	| FILL_COORDS); pElInfo; pElInfo = traverse_next(stack, pElInfo)) {
		if (pom <= iter) {
			cout<<(int)10*(pom/OneTenth)<<"%|";
			pom += OneTenth;
			//			pom = (int)10*(iter/OneTenth);
			//			cout<<pom<<"%|";
		}
		flush(cout);
		Vol = el_det(pElInfo);
		RestrictedElement = true;
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			pAlberta->GetGlobDof[i](pElInfo->el, pAlberta->FEspace[i]->admin,
					pAlberta->Dofs[i]);
			if (RestrictedElement)
				for (int j = 0; j < pAlberta->NumOfBasFct[i]; j++) {
					pAlberta->RestrictedDofs[i][j]
							=(DOF)pAlberta->pExtendedVertexIndex[i]->vec[pAlberta->Dofs[i][j]];
					if (pAlberta->RestrictedDofs[i][j] == -1) {
						RestrictedElement = false;
						break;
					}
				}
		}

		ComputeForElInfo(pElInfo);

		for (int i = 0; i<pAlberta->DimOfProb; i++)
			for (int j = 0; j<pAlberta->DimOfProb; j++) {
				if (pMatrixMask[i * pAlberta->DimOfProb + j]) {
					add_element_matrix(pAlberta->Matrix[i][j], Vol,
							pAlberta->NumOfBasFct[i], pAlberta->NumOfBasFct[j],
							pAlberta->Dofs[i], pAlberta->Dofs[j],
							(const REAL**)pAlberta->LocalMatrix[i][j], nil);
					if (RestrictedElement)
						add_element_matrix(pAlberta->RestrictedMatrix[i][j],
								Vol, pAlberta->NumOfBasFct[i],
								pAlberta->NumOfBasFct[j],
								pAlberta->RestrictedDofs[i],
								pAlberta->RestrictedDofs[j],
								(const REAL**)pAlberta->LocalMatrix[i][j], nil);
				}
			}
		iter++;
	}
	pAlberta->AssembleProblemMatrix(pMatrixMask,
			(const DOF_MATRIX ***)pAlberta->Matrix);
	//print_dof_matrix(pAlberta->ProblemMatrix);
	//WAIT_REALLY;
	if (pom <= iter)
		cout<<"100%|";
	cout<<endl;
	if (ActualLevel == NumOfLevels-1) // the highest level only
		pAlberta->PostProcessMatrix(1.);
	return 0;
}

/**
 Function calls user provided build on lowest level
 */
void CLevel::UserBuild() {
	if (User.pBuildDynamic) {
		User.pBuildDynamic(pAlberta->FEspace, pAlberta->Matrix);
		pAlberta->SumMatrices(pMatrixMask); // sum of static and dynamic parts
	}
	if (!User.pBuildDynamic && User.pBuildStatic) // second should be obvious
		pAlberta->AssembleProblemMatrix(pMatrixMask,
				(const DOF_MATRIX ***)pAlberta->StaticMatrix);
	else
		pAlberta->AssembleProblemMatrix(pMatrixMask,
				(const DOF_MATRIX ***)pAlberta->Matrix);
}

/**
 Same as Build. Only elements inside of domain are taken in account. Assemblage
 of system is not called.
 */
void CLevel::RestrictedUserBuild() {
	if (User.pBuildDynamic) {
		User.pBuildDynamic(pAlberta->RestrictedFEspace,
				pAlberta->RestrictedMatrix);
		pAlberta->SumRestrictedMatrices(pMatrixMask); // sum of static and dynamic parts
	}
}

/**
 Set up righthand side using function provided by user
 */
int CLevel::SetUpRHS() {

#ifdef USING_MPI
	if (ActualLevel==NumOfLevels-1) {
		if (MPIInfo.Rank==0)
			SetUpRHSMaster();
		else
			SetUpRHSSlave();
		return 0;
	}
#endif

	TRAVERSE_STACK *stack;
	int NumOfLeafEl = 0;
	int iter = 0;
	double pom = 0;

	NumOfLeafEl = GetNumOfLeafEl();
	double OneTenth = NumOfLeafEl/10.;

	for (int i = 0; i < pAlberta->DimOfProb; i++)
		dof_set(0.0, pAlberta->Fh[i]);

	stack = get_traverse_stack();
	pom = OneTenth;
	for (pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL
	| FILL_COORDS | FILL_BOUND); pElInfo; pElInfo = traverse_next(stack,
			pElInfo)) {
		if (pom <= iter) {
			cout<<(int)10*(pom/OneTenth)<<"%|";
			pom += OneTenth;
		}
		flush(cout);

		for (int i = 0; i < pAlberta->DimOfProb; i++)
			pAlberta->GetGlobDof[i](pElInfo->el, pAlberta->FEspace[i]->admin,
					pAlberta->Dofs[i]);

		SetUpRHSForOneElement(pElInfo);

		for (int i = 0; i<pAlberta->DimOfProb; i++)
			for (int j = 0; j<pAlberta->NumOfBasFct[i]; j++) {
				pAlberta->Fh[i]->vec[pAlberta->Dofs[i][j]]
						+= pAlberta->LocalFh[i][j];
			}

		iter++;
	}
	if (pom <= iter)
		cout<<"100%|";
	cout<<endl;

	return 0;
}

/**
 Problem solving on highest level. Other levels are assembled.
 */
void CLevel::Solve() {
	FUNCNAME("CLevel::Solve");

	INFO(9,9,"Solve started\n");

	SetUpRHS();
	DirichletBound();

	int idx=0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->ProblemRHS->vec[idx]=pAlberta->Fh[i]->vec[j];
			idx++;
		}
	}

	print_dof_matrix(pAlberta->ProblemMatrix);
	print_dof_real_vec(pAlberta->ProblemRHS);
	//WAIT_REALLY;
	//	User.pBuild(pAlberta->FEspace, pAlberta->Matrix);

	dof_set(0.0, pAlberta->ProblemSol);
	oem_solve_s(pAlberta->ProblemMatrix, pAlberta->ProblemRHS,
			pAlberta->ProblemSol, Solver.solver, Solver.tol, Solver.icon,
			Solver.restart, Solver.miter, Solver.info);
	//for (int j = 0; j < pAlberta->ProblemFEspace->admin->size_used; j++) cout<<pAlberta->ProblemSol->vec[j]<<"|";
	//WAIT_REALLY;
}

/**
 Dirichlet boundary condition used while interpolating solution from higher
 levels to the lower ones.
 */
REAL CLevel::Dirichlet(const REAL_D x/** world coordinates*/) {
	REAL Lambda[N_LAMBDA] = { 0 };
	REAL Val;
	world_to_coord(Aux.PrevElInfo, x, Lambda);
	Val
			= pPrevLevelStatic->pAlberta->FEspace[Aux.IndexOfFEspace]->bas_fcts->phi[Aux.IndexOfBasFct](Lambda);
	return Val;
}

void CLevel::SolveInt() {
	oem_solve_s(pAlberta->ProblemMatrix, pAlberta->ProblemRHS,
			pAlberta->ProblemUh[Aux.IndexOfFEspace][Aux.IndexOfBasFct],
			Solver.solver, Solver.tol, Solver.icon, Solver.restart,
			Solver.miter, Solver.info);
}

/**
 Interpolation of initial conditions going from higher level to lower.
 */
void CLevel::InterpolateFromPreviousLevel(const EL_INFO* pElInfo,
		int IndexOfFEspace, int IndexOfBF) {

	pPrevLevelStatic = pPrevLevel;
	Aux.PrevElInfo = pElInfo;
	Aux.IndexOfBasFct = IndexOfBF;
	Aux.IndexOfFEspace = IndexOfFEspace;

	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		dof_set(0.0, pAlberta->Fh[i]);
		dof_set(0.0, pAlberta->Uh[i]);//testuj vykonnost!!!
		dirichlet_bound((DIRICHLET_FUNC)&CLevel::Dirichlet, pAlberta->Fh[i],
				pAlberta->Uh[i], nil);
	}
	dof_set(0.0, pAlberta->ProblemRHS);
	dof_set(0.0, pAlberta->ProblemUh[Aux.IndexOfFEspace][Aux.IndexOfBasFct]);

	int idx=0;

	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->ProblemRHS->vec[idx]=pAlberta->Fh[i]->vec[j];
			pAlberta->ProblemUh[Aux.IndexOfFEspace][Aux.IndexOfBasFct]->vec[idx]
					=pAlberta->Uh[i]->vec[j];
			idx++;
		}
	}
	//for (int j = 0; j < pAlberta->ProblemFEspace->admin->size_used; j++) cout<<pAlberta->ProblemRHS->vec[j]<<"|";
	//print_dof_real_vec(pAlberta->ProblemRHS);
	//WAIT_REALLY;
}

/**
 Integration of 2 basis functions is same as xAy
 */
REAL CLevel::Integrate(int i, int j, int k, int l) {
	REAL ScaleFac, val;
	ScaleFac = 1./pow(ScaleFactor/2.0, pAlberta->DimOfWorld);
	dof_set(0.0, pAlberta->ProblemRHS); //pAlberta->Fh is auxiliary vector here
	dof_mv(NoTranspose, pAlberta->ProblemMatrix, pAlberta->ProblemUh[i][k],
			pAlberta->ProblemRHS);
	val = dof_dot(pAlberta->ProblemRHS, pAlberta->ProblemUh[j][l]);
	return val*ScaleFac;
}

/**
 Function returning solution
 */
DOF_REAL_VEC** CLevel::GetSolution() {
	int idx = 0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++) {
			pAlberta->Uh[i]->vec[j] = pAlberta->ProblemSol->vec[idx];
			idx++;
		}
	}
	return pAlberta->Uh;
}

/**
 Get number of leaf elements
 */
int CLevel::GetNumOfLeafEl() {
	TRAVERSE_STACK *stack = get_traverse_stack();
	int NumOfLeafEl = 0;
	for (pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL); pElInfo; pElInfo
			= traverse_next(stack, pElInfo))
		NumOfLeafEl++;
	// NumOfLeafEl == pAlberta->pMesh->n_elements
	free_traverse_stack(stack);
	return NumOfLeafEl;
}

void CLevel::ComputeForElInfo(const EL_INFO * pElInfo) {
	REAL_D x;

	pAlberta->ClearLocalMatrix();

	//parametric meshes not supported!
	for (int iq = 0; iq < pAlberta->pQuad->n_points; iq++) {
		coord_to_world(pElInfo, pAlberta->pQuad->lambda[iq], x);
		pNextLevel->pAlberta->TransformMeshes((REAL *)x);
		if (ActualLevel>1) {
			pNextLevel->Assemble();
		} else {
			pNextLevel->pAlberta->ClearMatrix();
			pNextLevel->pAlberta->ClearRestrictedMatrix();
			pNextLevel->UserBuild();
			pNextLevel->RestrictedUserBuild();
		} // at this point we surely have DOF_MATRIX assembled at a lower level
		for (int j = 0; j<pAlberta->DimOfProb; j++) {
			for (int k = 0; k<pAlberta->NumOfBasFct[j]; k++) {
				pNextLevel->InterpolateFromPreviousLevel(pElInfo, j, k);
				pNextLevel->SolveInt();
				if (ActualLevel == 1)
					pNextLevel->PostprocessMicroSolution();
			}
		}
		for (int i = 0; i<pAlberta->DimOfProb; i++)
			for (int j = 0; j<pAlberta->DimOfProb; j++) {
				if (pMatrixMask[i * pAlberta->DimOfProb + j]) {
					pNextLevel->ExtendRestrictedMatrix(i, j);
					for (int k = 0; k<pAlberta->NumOfBasFct[i]; k++) {
						for (int l = 0; l<pAlberta->NumOfBasFct[j]; l++) {
							pAlberta->LocalMatrix[i][j][k][l]
									+= pAlberta->pQuad->w[iq]
											*pNextLevel->Integrate(i, j, k, l);//optimize
						}
					}
				}
			}
	}
}

void CLevel::ExtendRestrictedMatrix(int i, int j) {
	if (ActualLevel == 0 && !User.pBuildDynamic)
		pAlberta->ExtendRestrictedMatrix(i, j,
				(const DOF_MATRIX ***)pAlberta->StaticRestrictedMatrix);
	else
		pAlberta->ExtendRestrictedMatrix(i, j,
				(const DOF_MATRIX ***)pAlberta->RestrictedMatrix);
}

void CLevel::SetUpScaleFactor(REAL Prev) {
	FUNCNAME("CLevel::SetUpScaleFactor");

	REAL Radius, MinRadius;
	REAL_D dx;
	TRAVERSE_STACK *stack;
	const EL_INFO* pElInfo;

	ScaleFactor *= Prev;
	for (int i = 0; i < pAlberta->DimOfWorld; i++)
		// transformation is only scaling not transition
		dx[i] = ScaleFactor/2.0;

	pAlberta->ScaleFactor = ScaleFactor;
	pAlberta->TransformMeshes(dx); // on the highest level is transformation permanent


	stack = get_traverse_stack();

	pElInfo = traverse_first(stack, pAlberta->pMesh, -1, CALL_LEAF_EL
	| FILL_COORDS);
	MinRadius = pAlberta->GetInnerRadius(pElInfo);

	for (pElInfo = traverse_next(stack, pElInfo); pElInfo; pElInfo
			= traverse_next(stack, pElInfo)) {
		Radius = pAlberta->GetInnerRadius(pElInfo);
		if (MinRadius > Radius)
			MinRadius = Radius;
	}
	free_traverse_stack(stack);

	if (pNextLevel)
		pNextLevel->SetUpScaleFactor(MinRadius);

#ifdef USING_MPI
	if (MPIInfo.Rank != 0) // status info only on master
		return;
#endif

	INFO(1,1,"ActualLevel = %d       MeshDiameter = %f\n",ActualLevel,MinRadius);
}

void CLevel::SetUpRHSForOneElement(const EL_INFO *pElInfo) {
	for (int i = 0; i < pAlberta->DimOfProb; i++)
		for (int j = 0; j < pAlberta->NumOfBasFct[i]; j++)
			pAlberta->LocalFh[i][j] = 0;
	UserSetUpRHS(pElInfo);
}

void CLevel::UserSetUpRHS(const EL_INFO *pElInfo) {
	User.pSetUpRHS(pElInfo, pAlberta->FEspace, pAlberta->LocalFh);
}

REAL CLevel::UserDirichletBound(const REAL_D x) {
	return pDirBoundCond(x, SpaceIndex);
}

void CLevel::DirichletBound() {
	if (pDirBoundCond)
		for (int i = 0; i < pAlberta->DimOfProb; i++) {
			SpaceIndex = i;
			dirichlet_bound(&CLevel::UserDirichletBound, pAlberta->Fh[i],
					pAlberta->Uh[i], nil);
		}
}

void CLevel::PrintSolution(int i) {
	for (int j = 0; j < pAlberta->FEspace[i]->admin->size_used; j++)
		cout<<pAlberta->Uh[i]->vec[j]<<"|";
	cout<<endl;
}

void CLevel::PostprocessMicroSolution() {
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
	User.pMicroSolPostprocess(pAlberta->FEspace, pAlberta->Uh);
}

void CLevel::UserBuildStatic() {
	REAL_D x;

	if (!User.pBuildStatic)
		return;
	if (!pAlberta->StaticMatrix)
		pAlberta->CreateStaticMatrices();

	for (int i = 0; i < pAlberta->DimOfWorld; i++)
		x[i] = (pAlberta->x[i] + pAlberta->y[i])/2.;
	pAlberta->TransformMeshes(x);
	User.pBuildStatic(pAlberta->FEspace, pAlberta->StaticMatrix);
	User.pBuildStatic(pAlberta->RestrictedFEspace,
			pAlberta->StaticRestrictedMatrix);
}

void CLevel::Solve(bool WantAssemble) {
	FUNCNAME("CLevel::Solve(bool WantAssemble)");

	if (!IsAssembled || WantAssemble) { // I want to assemble
#ifdef USING_MPI
		if (MPIInfo.Rank == 0)
#endif
			if (!IsAssembled)
				INFO(9,9,"First assembling of system matrix\n");

		Assemble();
		IsAssembled = true;
	}

	Solve();
	L2Error = CalcL2Error();
}

void CLevel::SaveSolutionToVTK(const char *pFileName) {
	pSolution = GetSolution();
	pSaveVTK->SaveVTKfile((const DOF_REAL_VEC **)pSolution, pFileName);
}

REAL CLevel::UserLocalExact(const REAL_D x) {
	return pExact(x, ActComponentIdx);
}

REAL CLevel::LocalExact(const REAL_D x) {
	return UserLocalExact(x);
}

REAL CLevel::CalcL2Error() {
	FUNCNAME("CLevel::CalcL2Error()");

	REAL RetVal, Err;
	pSolution = GetSolution();
	RetVal = 0;
	for (int i = 0; i < pAlberta->DimOfProb; i++) {
		ActComponentIdx = i;
		Err
				= L2_err(&CLevel::UserLocalExact, pSolution[i], nil, true, nil,
						nil);
		RetVal += Err*Err;
	}

#ifdef USING_MPI
	if (MPIInfo.Rank == 0)
#endif
		INFO(9,9,"L2error local = %f\n",pow(RetVal, 0.5));
	return pow(RetVal, 0.5);
}

