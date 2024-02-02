#include "CAlbertaInterface.h"

#define VERT_IND(dim,i,j) ((i)*N_VERTICES(dim)+(j))
#define NEIGH_IND(dim,i,j) ((i)*N_NEIGH(dim)+(j))

using namespace std;

void CAlbertaInterface::CleanupWriteMacro(MACRO_DATA *data,
		DOF_INT_VEC *dof_vert_ind, TRAVERSE_STACK *stack) {

	if (data)
		free_macro_data(data);
	free_dof_int_vec(dof_vert_ind);
	free_traverse_stack(stack);
	return;
}

int CAlbertaInterface::IsInside(const EL_INFO* el_info) {
	/* returns 0 if element is outside internal square, 1 otherwise */
	for (int i = 0; i < DimOfWorld+1; i++)
		for (int j = 0; j < DimOfWorld; j++)
			if ((el_info->coord[i][j] < x[j]) || (y[j] < el_info->coord[i][j]))
				return 0;
	return 1;
}

int CAlbertaInterface::MeshToMacroData(MESH *pMesh) {
	TRAVERSE_STACK *stack;
	FLAGS fill_flag= CALL_LEAF_EL|FILL_COORDS|FILL_BOUND|FILL_NEIGH;
	PARAMETRIC *parametric;
	const DOF_ADMIN *admin;
	FE_SPACE fe_space = { "write fe_space", nil, nil };
	const EL_INFO *el_info;
	DOF_INT_VEC *dof_vert_ind;
	int dim = pMesh->dim, n0, ne, nv, i, j, *vert_ind= nil;
	U_CHAR write_el_type;
	static const REAL vertex_bary[N_VERTICES_MAX][N_LAMBDA] = { { 1.0, 0.0,
			0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 }, { 0.0,
			0.0, 0.0, 1.0 } };

	admin = get_vertex_admin(pMesh);

	n0 = admin->n0_dof[VERTEX];
	fe_space.admin = admin;

	parametric = pMesh->parametric;

	dof_vert_ind = get_dof_int_vec("vertex indices", &fe_space);
	GET_DOF_VEC(vert_ind, dof_vert_ind);
	FOR_ALL_DOFS(admin, vert_ind[dof] = -1);

	nv = ne = 0;
	write_el_type = false;

	stack = get_traverse_stack();
	///< this cycle is to compute the number of inner vertices and elements  
	for (el_info = traverse_first(stack, pMesh, -1, CALL_LEAF_EL | FILL_COORDS); el_info; el_info
			= traverse_next(stack, el_info)) {

		if (!IsInside(el_info))
			continue;

		if (parametric) {
			parametric->init_element(el_info, parametric);
			parametric->coord_to_world(el_info, nil, N_VERTICES(dim), vertex_bary, (REAL_D *)el_info->coord);
		}

		for (i = 0; i < N_VERTICES(dim); i++) {
			if (vert_ind[el_info->el->dof[i][n0]] == -1) {
				vert_ind[el_info->el->dof[i][n0]] = 1;
				nv++;
			}
		}
		ne++;
	}

	if (pData)
		free_macro_data(pData);
	pData = alloc_macro_data(dim, nv, ne, 0);
	nv = ne = 0;
	FOR_ALL_DOFS(admin, vert_ind[dof] = -1);

	/*--------------------------------------------------------------------------*/
	/* The first pass counts elements and vertices, checks these against the    */
	/* entries of pMesh->n_elements, pMesh->n_vertices, and fills pData->coords.   */
	/* A check on whether an element has nonzero el_type is also done.          */
	/*--------------------------------------------------------------------------*/

	for (el_info = traverse_first(stack, pMesh, -1, CALL_LEAF_EL | FILL_COORDS); el_info; el_info
			= traverse_next(stack, el_info)) {

		if (!IsInside(el_info))
			continue;

		if (parametric) {
			parametric->init_element(el_info, parametric);
			parametric->coord_to_world(el_info, nil, N_VERTICES(dim), vertex_bary, (REAL_D *)el_info->coord);
		}

		for (i = 0; i < N_VERTICES(dim); i++) {
			if (vert_ind[el_info->el->dof[i][n0]] == -1) {
				/*--------------------------------------------------------------------------*/
				/* assign a global index to each vertex                                     */
				/*--------------------------------------------------------------------------*/
				vert_ind[el_info->el->dof[i][n0]] = nv;

				for (j = 0; j < DIM_OF_WORLD; j++)
					pData->coords[nv][j] = el_info->coord[i][j];

				nv++;

				if (nv > pMesh->n_vertices) {
					CleanupWriteMacro(pData, dof_vert_ind, stack);
					ERROR("pMesh %s: n_vertices (==%d) is too small! Writing aborted\n",
							pMesh->name, pMesh->n_vertices);
					return -1;
				}
			}
		}

		ne++;

		if (ne > pMesh->n_elements) {
			CleanupWriteMacro(pData, dof_vert_ind, stack);
			ERROR("pMesh %s: n_elements (==%d) is too small! Writing aborted\n",
					pMesh->name, pMesh->n_elements);
			return -1;
		}
		if (dim == 3 && el_info->el_type)
			write_el_type = true;
	}

	if (dim > 0)
		pData->boundary = MEM_ALLOC(ne*N_NEIGH(dim), S_CHAR);

	if (write_el_type)
		pData->el_type = MEM_ALLOC(ne, U_CHAR);

	ne = 0;

	/*--------------------------------------------------------------------------*/
	/* The second pass assigns mel_vertices, boundary, and if necessary el_type */
	/*--------------------------------------------------------------------------*/
	for (el_info = traverse_first(stack, pMesh, -1, fill_flag); el_info; el_info
			= traverse_next(stack, el_info)) {

		if (!IsInside(el_info))
			continue;

		for (i = 0; i < N_VERTICES(dim); i++)
			pData->mel_vertices[VERT_IND(dim,ne,i)] = vert_ind[el_info->el->dof[i][n0]];

		if (dim > 0)
			for (i = 0; i < N_NEIGH(dim); i++)
				switch (dim) {
				case 1:
					pData->boundary[NEIGH_IND(dim,ne,i)] = el_info->vertex_bound[1-i];
					break;
				case 2:
					pData->boundary[NEIGH_IND(dim,ne,i)] = el_info->edge_bound[i];
					break;
				case 3:
					pData->boundary[NEIGH_IND(dim,ne,i)] = el_info->face_bound[i];
				}

		if (write_el_type)
			pData->el_type[ne] = el_info->el_type;

		++ne;
	}

	/*--------------------------------------------------------------------------*/
	/* Finally, we compute neighbour information. This seems to be the easiest  */
	/* solution, since neighbor information in ALBERTA is only available as     */
	/* pointers.                                                                */
	/*--------------------------------------------------------------------------*/

	if (dim > 0)
		compute_neigh_fast(pData);
	CleanupWriteMacro(nil, dof_vert_ind, stack);
	return 0;
}

MESH* CAlbertaInterface::GetRestrictedMesh(REAL *xp, REAL *yp, MESH* pMesh) {
	DimOfWorld=pMesh->dim;
	for (int i=0; i<DimOfWorld; i++) {
		x[i]=xp[i];
		y[i]=yp[i];
	}
	MeshToMacroData(pMesh);
	return GetInteriorMesh();
}

#include "macro_1d.c"
#if DIM_OF_WORLD > 1
#include "macro_2d.c"
#if DIM_OF_WORLD > 2
#include "macro_3d.c"
#endif
#endif

MESH* CAlbertaInterface::GetInteriorMesh() {
	MESH* pMesh;
	int i, j;
	for (i=0; i<pData->n_macro_elements; i++)
		/* pretend it is Dirichlet to fool GET_MESH routine */
		for (j=0; j<N_NEIGH(DimOfWorld); j++)
			if (pData->neigh[i*N_NEIGH(DimOfWorld)+j]==-1)
				pData->boundary[i*N_NEIGH(DimOfWorld)+j]=DIRICHLET;
	pMesh=GET_MESH(DimOfWorld, "RestrictedMESH", pData, nil);
	for (i=0; i<pData->n_macro_elements; i++)
		/* all BOUNDARY is interior */
		for (j=0; j<N_NEIGH(DimOfWorld); j++)
			if (pData->neigh[i*N_NEIGH(DimOfWorld)+j]==-1)
				pData->boundary[i*N_NEIGH(DimOfWorld)+j]=INTERIOR;
	switch (DimOfWorld) {
	case 1:
		fill_bound_info_1d(pMesh, pData);
		break;
#if DIM_OF_WORLD > 1
		case 2:
		fill_bound_info_2d(pMesh, pData);
		break;
#if DIM_OF_WORLD > 2
		case 3:
		fill_bound_info_3d(pMesh, pData);
		fill_more_bound_info_3d(pMesh, pData);
		break;
#endif
#endif
	default:
		ERROR_EXIT("Illegal dimension %d!\n", DimOfWorld);
	}
	return pMesh;
}

MACRO_DATA * CAlbertaInterface::GetReferenceMacroData(int DimOfWorld) {
	MACRO_DATA *macro_data;
	int nv, ne;
	ne=nv=0;
	switch (DimOfWorld) {
	case 1:
		ne=1;
		nv=2;
		break;
	case 2:
		ne=4;
		nv=5;
		break;
	case 3:
		ne=6;
		nv=8;
		break;
	}
	if (!(ne*nv))
		return NULL;
	macro_data = alloc_macro_data(DimOfWorld, nv, ne, 0);

	// element coordinates
	switch (DimOfWorld) {
	case 1:
		macro_data->coords[0][0]=0.0;
		macro_data->coords[1][0]=1.0;
		break;
	case 2:
		macro_data->coords[0][0]=0.0;
		macro_data->coords[0][1]=0.0;
		macro_data->coords[1][0]=1.0;
		macro_data->coords[1][1]=0.0;
		macro_data->coords[2][0]=1.0;
		macro_data->coords[2][1]=1.0;
		macro_data->coords[3][0]=0.0;
		macro_data->coords[3][1]=1.0;
		macro_data->coords[4][0]=0.5;
		macro_data->coords[4][1]=0.5;
		break;
	case 3:
		macro_data->coords[0][0]=0.0;
		macro_data->coords[0][1]=0.0;
		macro_data->coords[0][2]=0.0;
		macro_data->coords[1][0]=1.0;
		macro_data->coords[1][1]=0.0;
		macro_data->coords[1][2]=0.0;
		macro_data->coords[2][0]=0.0;
		macro_data->coords[2][1]=0.0;
		macro_data->coords[2][2]=1.0;
		macro_data->coords[3][0]=1.0;
		macro_data->coords[3][1]=0.0;
		macro_data->coords[3][2]=1.0;
		macro_data->coords[4][0]=1.0;
		macro_data->coords[4][1]=1.0;
		macro_data->coords[4][2]=0.0;
		macro_data->coords[5][0]=1.0;
		macro_data->coords[5][1]=1.0;
		macro_data->coords[5][2]=1.0;
		macro_data->coords[6][0]=0.0;
		macro_data->coords[6][1]=1.0;
		macro_data->coords[6][2]=0.0;
		macro_data->coords[7][0]=0.0;
		macro_data->coords[7][1]=1.0;
		macro_data->coords[7][2]=1.0;
		break;
	}

	// macro element vertices
	switch (DimOfWorld) {
	case 1:
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+1]=1;
		break;
	case 2:
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+1]=1;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+2]=4;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+0]=1;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+1]=2;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+2]=4;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+0]=2;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+1]=3;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+2]=4;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+0]=3;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+1]=0;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+2]=4;
		break;
	case 3:
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+2]=4;
		macro_data->mel_vertices[0*N_VERTICES(DimOfWorld)+3]=1;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+2]=3;
		macro_data->mel_vertices[1*N_VERTICES(DimOfWorld)+3]=1;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+2]=3;
		macro_data->mel_vertices[2*N_VERTICES(DimOfWorld)+3]=2;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+2]=4;
		macro_data->mel_vertices[3*N_VERTICES(DimOfWorld)+3]=6;
		macro_data->mel_vertices[4*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[4*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[4*N_VERTICES(DimOfWorld)+2]=7;
		macro_data->mel_vertices[4*N_VERTICES(DimOfWorld)+3]=6;
		macro_data->mel_vertices[5*N_VERTICES(DimOfWorld)+0]=0;
		macro_data->mel_vertices[5*N_VERTICES(DimOfWorld)+1]=5;
		macro_data->mel_vertices[5*N_VERTICES(DimOfWorld)+2]=7;
		macro_data->mel_vertices[5*N_VERTICES(DimOfWorld)+3]=2;
		break;
	}

	// boundary
	macro_data->boundary = MEM_ALLOC(ne*N_NEIGH(DimOfWorld), S_CHAR);
	switch (DimOfWorld) {
	case 1:
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+1]=1;
		break;
	case 2:
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+0]=0;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+1]=0;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+2]=1;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+0]=0;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+1]=0;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+2]=1;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+0]=0;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+1]=0;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+2]=1;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+0]=0;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+1]=0;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+2]=1;
		break;
	case 3:
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[0*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[1*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[2*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[3*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->boundary[4*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[4*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[4*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[4*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->boundary[5*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->boundary[5*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->boundary[5*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->boundary[5*N_NEIGH(DimOfWorld)+3]=0;
		break;
	}

	macro_data->neigh = MEM_ALLOC(ne*N_NEIGH(DimOfWorld), int);
	switch (DimOfWorld) {
	case 1:
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+1]=-1;
		break;
	case 2:
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+0]=1;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+1]=3;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+2]=-1;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+0]=2;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+1]=0;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+2]=-1;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+0]=3;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+1]=1;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+2]=-1;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+0]=0;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+1]=2;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+2]=-1;
		break;
	case 3:
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+2]=1;
		macro_data->neigh[0*N_NEIGH(DimOfWorld)+3]=3;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+2]=0;
		macro_data->neigh[1*N_NEIGH(DimOfWorld)+3]=2;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+2]=5;
		macro_data->neigh[2*N_NEIGH(DimOfWorld)+3]=1;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+2]=4;
		macro_data->neigh[3*N_NEIGH(DimOfWorld)+3]=0;
		macro_data->neigh[4*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[4*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[4*N_NEIGH(DimOfWorld)+2]=3;
		macro_data->neigh[4*N_NEIGH(DimOfWorld)+3]=5;
		macro_data->neigh[5*N_NEIGH(DimOfWorld)+0]=-1;
		macro_data->neigh[5*N_NEIGH(DimOfWorld)+1]=-1;
		macro_data->neigh[5*N_NEIGH(DimOfWorld)+2]=2;
		macro_data->neigh[5*N_NEIGH(DimOfWorld)+3]=4;
		break;
	}
	return (macro_data);
}

int CAlbertaInterface::GetProblemDof(int NumOfSpace, int Dof) {
	return IdxOfDofs[NumOfSpace]+Dof;
}

void CAlbertaInterface::Initialize(
		CUserData * Data/** Data provided by user */, int ActLevel, int Rank) {
	ActualLevel = ActLevel;
	NumOfLevels = Data->NumOfLevels;
	ScaleFactor = Data->pScaleFactor[NumOfLevels-ActualLevel-1];
	//Alberta
	DimOfWorld = Data->DimOfWorld;
	DimOfProb = Data->DimOfProb;
	NumOfRefin = Data->pNumOfRefin[NumOfLevels-ActualLevel-1];
	QuadDeg = 1;
	pQuad = get_quadrature(DimOfWorld, QuadDeg);
	LagDeg = (int *)malloc(DimOfProb*sizeof(int));
	BoundCond = Data->BoundCond;

	StaticMatrix = NULL;
	StaticRestrictedMatrix = NULL;

	if (ActualLevel>0) {
		for (int i = 0; i < DimOfProb; i++)
			LagDeg[i] = 1;
	} else {
		for (int i = 0; i < DimOfProb; i++)
			LagDeg[i] = Data->pFEdeg[i];
	}

	InitGeometry(Data->pMacroFile, Rank);

	pBoundaryEl = NULL;
	LinearFEspace = NULL;
	ConstrMatrix = NULL;
	ConstrMatrixT = NULL;

	pAllBoundaryEl = NULL;
	pOppositeDOFs = NULL;
	if (BoundCond==1) { //if we have the periodic boundary condition
		GetNumOfBoundEl();
		pBoundaryEl=(MACRO_EL****)malloc(DimOfWorld*sizeof(MACRO_EL***));
		for (int i = 0; i < DimOfWorld; i++) {
			pBoundaryEl[i] = (MACRO_EL***)malloc(2*sizeof(MACRO_EL**));
			pBoundaryEl[i][0] = (MACRO_EL**)malloc(NumOfBoundEl
					*sizeof(MACRO_EL*));
			pBoundaryEl[i][1] = (MACRO_EL**)malloc(NumOfBoundEl
					*sizeof(MACRO_EL*));
		}
		pAllBoundaryEl = (MACRO_EL**)malloc(NumOfBoundEl*sizeof(MACRO_EL*));

		const BAS_FCTS *Lag = get_lagrange(DimOfWorld, 1);
		LinearFEspace = get_fe_space(pMesh, Lag->name, nil, Lag, 
		false);
		// LinearFEspace and pOppositeDOFs are needed for this 2 functions to work
		pOppositeDOFs=get_dof_int_vec("OppositeDOFs", LinearFEspace);
		memset(pOppositeDOFs->vec, -1, pOppositeDOFs->size*sizeof(int));

		FillPrimaryBoundElArray();
		FillOppositeBoundElArray();

		ConstrMatrix = (DOF_MATRIX**)malloc(DimOfProb*sizeof(DOF_MATRIX*));
		ConstrMatrixT = (DOF_MATRIX**)malloc(DimOfProb*sizeof(DOF_MATRIX*));
		for (int i = 0; i < DimOfProb; i++) {
			ConstrMatrix[i] = get_dof_matrix("D", LinearFEspace, FEspace[i]);
			ConstrMatrixT[i] = get_dof_matrix("DT", FEspace[i], LinearFEspace);
		}
	}

	NumOfBasFct = (int*)malloc(DimOfProb*sizeof(int));
	GetGlobDof
			= (const DOF *(**)(const EL *, const DOF_ADMIN *, DOF *))malloc(DimOfProb
					*sizeof(const DOF *(*)(const EL *, const DOF_ADMIN *, DOF *)));
	GetBound = (const S_CHAR *(**)(const EL_INFO *, S_CHAR *))malloc(DimOfProb
			*sizeof(const S_CHAR *(*)(const EL_INFO *, S_CHAR *)));
	Dofs = (DOF**)malloc(DimOfProb*sizeof(DOF*));
	RestrictedDofs = (DOF**)malloc(DimOfProb*sizeof(DOF*));
	for (int i = 0; i < DimOfProb; i++) {
		NumOfBasFct[i] = FEspace[i]->bas_fcts->n_bas_fcts;
		GetGlobDof[i] = FEspace[i]->bas_fcts->get_dof_indices;
		GetBound[i] = FEspace[i]->bas_fcts->get_bound;
		Dofs[i] = (DOF*)malloc(NumOfBasFct[i]*sizeof(DOF));
		RestrictedDofs[i] = (DOF*)malloc(NumOfBasFct[i]*sizeof(DOF));
	}
	Matrix = (DOF_MATRIX***)malloc(DimOfProb *sizeof(DOF_MATRIX**));
	for (int i = 0; i < DimOfProb; i++) {
		Matrix[i] = (DOF_MATRIX**)malloc(DimOfProb *sizeof(DOF_MATRIX*));
		for (int j = 0; j < DimOfProb; j++) {
			Matrix[i][j] = get_dof_matrix("A", FEspace[i], FEspace[j]);
		}
	}
	RestrictedMatrix = (DOF_MATRIX***)malloc(DimOfProb *sizeof(DOF_MATRIX**));
	for (int i = 0; i < DimOfProb; i++) {
		RestrictedMatrix[i] = (DOF_MATRIX**)malloc(DimOfProb
				*sizeof(DOF_MATRIX*));
		for (int j = 0; j < DimOfProb; j++) {
			RestrictedMatrix[i][j] = get_dof_matrix("ARes",
					RestrictedFEspace[i], RestrictedFEspace[j]);
		}
	}

	pRestrictedVertexIndex = (DOF_REAL_VEC**)malloc(DimOfProb
			*sizeof(DOF_REAL_VEC*));
	for (int i = 0; i < DimOfProb; i++)
		pRestrictedVertexIndex[i] = get_dof_real_vec("RestrictedVectorMapping",
				RestrictedFEspace[i]);

	pExtendedVertexIndex = (DOF_REAL_VEC**)malloc(DimOfProb
			*sizeof(DOF_REAL_VEC*));
	for (int i = 0; i < DimOfProb; i++)
		pExtendedVertexIndex[i] = get_dof_real_vec("ExtendedVectorMapping",
				FEspace[i]);

	LocalMatrix=NULL;
	if (ActualLevel>0) {
		LocalMatrix = (REAL****)malloc(DimOfProb*sizeof(REAL***));
		for (int i = 0; i < DimOfProb; i++) {
			LocalMatrix[i] = (REAL***)malloc(DimOfProb*sizeof(REAL**));
			for (int j = 0; j < DimOfProb; j++) {
				LocalMatrix[i][j] = (REAL**)malloc(NumOfBasFct[i]
						*sizeof(REAL*));
				for (int k = 0; k < NumOfBasFct[i]; k++) {
					LocalMatrix[i][j][k] = (REAL*)malloc(NumOfBasFct[j]
							*sizeof(REAL));
				}
			}
		}
	}

	Uh = (DOF_REAL_VEC**)malloc(DimOfProb *sizeof(DOF_REAL_VEC*));
	for (int i = 0; i < DimOfProb; i++)
		Uh[i] = get_dof_real_vec("Uh", FEspace[i]);

	LocalFh = NULL;
	if (ActualLevel==NumOfLevels-1) {
		LocalFh = (REAL**)malloc(DimOfProb*sizeof(REAL*));
		for (int i = 0; i < DimOfProb; i++)
			LocalFh[i] = (REAL*)malloc(NumOfBasFct[i]*sizeof(REAL));
	}

	const int ndof[4]= { 0, 0, 0, 0 }; //HAHAHA!!!
	ProblemFEspace = get_fe_space(pMesh, "ProblemFEspace", ndof, nil, false);
	int FEsize = 0;
	for (int i = 0; i < DimOfProb; i++)
		FEsize += FEspace[i]->admin->size_used;

	DOF* pom;
	pom = (DOF*)&ProblemFEspace->admin->size_used;
	*pom = FEsize;
	pom = (DOF*)&ProblemFEspace->admin->size;
	*pom = FEsize;
	pom = (DOF*)&ProblemFEspace->admin->used_count;
	*pom = FEsize;

	ProblemMatrix = get_dof_matrix("ProblemMatrix", ProblemFEspace, nil);
	ProblemRHS = get_dof_real_vec("ProblemRHS", ProblemFEspace);
	ProblemSol = get_dof_real_vec("ProblemSol", ProblemFEspace);

	Fh = (DOF_REAL_VEC**)malloc(DimOfProb*sizeof(DOF_REAL_VEC*));
	for (int i = 0; i < DimOfProb; i++)
		Fh[i] = get_dof_real_vec("Fh", FEspace[i]);

	IdxOfDofs = new int[DimOfProb];
	IdxOfDofs[0] = 0;
	for (int i = 1; i < DimOfProb; i++)
		IdxOfDofs[i]=IdxOfDofs[i-1] +FEspace[i]->admin->size_used;
}

void CAlbertaInterface::CreateStaticMatrices() {
	StaticMatrix = (DOF_MATRIX***)malloc(DimOfProb *sizeof(DOF_MATRIX**));
	for (int i = 0; i < DimOfProb; i++) {
		StaticMatrix[i] = (DOF_MATRIX**)malloc(DimOfProb *sizeof(DOF_MATRIX*));
		for (int j = 0; j < DimOfProb; j++) {
			StaticMatrix[i][j] = get_dof_matrix("StaticA", FEspace[i],
					FEspace[j]);
		}
	}
	StaticRestrictedMatrix = (DOF_MATRIX***)malloc(DimOfProb
			*sizeof(DOF_MATRIX**));
	for (int i = 0; i < DimOfProb; i++) {
		StaticRestrictedMatrix[i] = (DOF_MATRIX**)malloc(DimOfProb
				*sizeof(DOF_MATRIX*));
		for (int j = 0; j < DimOfProb; j++) {
			StaticRestrictedMatrix[i][j] = get_dof_matrix("StaticARes",
					RestrictedFEspace[i], RestrictedFEspace[j]);
		}
	}
}

/**
 Function initializating geometry using user provided macro file (highest level only).
 All elements are macroelements.
 */
void CAlbertaInterface::InitGeometry(const char *MacroFile, int Rank) {
	MACRO_DATA* pData= NULL;
	char Name[100];
	sprintf(Name, "macro_tmp_file_%d.txt", Rank);
	if (ActualLevel == NumOfLevels-1)
		pData = read_macro(MacroFile);
	else
		pData = GetReferenceMacroData(DimOfWorld);
	pMesh = GET_MESH(DimOfWorld, "ALBERTA mesh", pData, nil);
	free_macro_data(pData);

	Lagrange = (const BAS_FCTS**)malloc(DimOfProb *sizeof(BAS_FCTS*));
	for (int i = 0; i < DimOfProb; i++) {
		Lagrange[i] = get_lagrange(DimOfWorld, LagDeg[i]);
	}
	get_fe_space(pMesh, Lagrange[0]->name, nil, Lagrange[0], 
	false);
	global_refine(pMesh, NumOfRefin*DimOfWorld);
	write_macro(pMesh, Name);
	fflush(NULL);
	free_mesh(pMesh);
	pData = read_macro(Name);
	pMesh = GET_MESH(DimOfWorld, "ALBERTA mesh", pData, nil);
	free_macro_data(pData);

	FEspace = (const FE_SPACE**)malloc(DimOfProb *sizeof(FE_SPACE*));
	for (int i = 0; i < DimOfProb; i++) {
		FEspace[i] = get_fe_space(pMesh, Lagrange[i]->name, nil, Lagrange[i], 
		false);
	}
	pVertices = new REAL[DimOfWorld*pMesh->n_elements*N_VERTICES(DimOfWorld)];
	for (int i = 0; i<pMesh->n_elements; i++)
		for (int j = 0; j < N_VERTICES(DimOfWorld); j++)
			for (int k = 0; k<DimOfWorld; k++)
				pVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k]
						= pMesh->macro_els[i].coord[j][k];
	unlink(Name);
	//now we initialize RestrictedFESpace...
	REAL_D x, y;
	for (int i = 0; i<DimOfWorld; i++) {
		x[i] = 0.25;
		y[i] = 0.75;
	}
	pRestrictedMesh = GetRestrictedMesh(x, y, pMesh);
	RestrictedFEspace = (const FE_SPACE**)malloc(DimOfProb *sizeof(FE_SPACE*));
	for (int i = 0; i < DimOfProb; i++) {
		RestrictedFEspace[i] = get_fe_space(pRestrictedMesh, Lagrange[i]->name, 
		nil, Lagrange[i], false);
	}
	pRestrictedVertices
			= new REAL[DimOfWorld*pRestrictedMesh->n_elements*N_VERTICES(DimOfWorld)];
	for (int i = 0; i<pRestrictedMesh->n_elements; i++)
		for (int j = 0; j < N_VERTICES(DimOfWorld); j++)
			for (int k = 0; k<DimOfWorld; k++)
				pRestrictedVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k]
						= pRestrictedMesh->macro_els[i].coord[j][k];
}

void CAlbertaInterface::Free() {
	if (pData) {
		free_macro_data(pData);
		pData = NULL;
	}
	if (pVertices) {
		delete []pVertices;
		pVertices = NULL;
	}
	if (pRestrictedVertices) {
		delete []pRestrictedVertices;
		pRestrictedVertices = NULL;
	}
	if (Dofs) {
		for (int i = 0; i < DimOfProb; i++) {
			free(Dofs[i]);
		}
		free(Dofs);
		Dofs = NULL;
	}
	if (RestrictedDofs) {
		for (int i = 0; i < DimOfProb; i++) {
			free(RestrictedDofs[i]);
		}
		free(RestrictedDofs);
		RestrictedDofs = NULL;
	}
	if (Lagrange) {
		free(Lagrange);
		Lagrange = NULL;
	}
	if (LagDeg) {
		free(LagDeg);
		LagDeg = NULL;
	}
	if (GetGlobDof) {
		free(GetGlobDof);
		GetGlobDof = NULL;
	}
	if (GetBound) {
		free(GetBound);
		GetBound = NULL;
	}
	if (LocalMatrix) {
		for (int k = 0; k < DimOfProb; k++) {
			for (int i = 0; i < DimOfProb; i++) {
				for (int j = 0; j < NumOfBasFct[i]; j++)
					free(LocalMatrix[k][i][j]);
				free(LocalMatrix[k][i]);
			}
			free(LocalMatrix[k]);
		}
		free(LocalMatrix);
		LocalMatrix = NULL;
	}
	if (LocalFh) {
		for (int i = 0; i < DimOfProb; i++)
			free(LocalFh[i]);
		free(LocalFh);
		LocalFh = NULL;
	}
	if (Matrix) {
		for (int i = 0; i < DimOfProb; i++) {
			for (int j = 0; j < DimOfProb; j++) {
				clear_dof_matrix(Matrix[i][j]);
			}
			free(Matrix[i]);
		}
		free(Matrix);
		Matrix = NULL;
	}
	if (RestrictedMatrix) {
		for (int i = 0; i < DimOfProb; i++) {
			for (int j = 0; j < DimOfProb; j++) {
				clear_dof_matrix(RestrictedMatrix[i][j]);
			}
			free(RestrictedMatrix[i]);
		}
		free(RestrictedMatrix);
		RestrictedMatrix = NULL;
	}
	if (StaticMatrix) {
		for (int i = 0; i < DimOfProb; i++) {
			for (int j = 0; j < DimOfProb; j++) {
				clear_dof_matrix(StaticMatrix[i][j]);
			}
			free(StaticMatrix[i]);
		}
		free(StaticMatrix);
		StaticMatrix = NULL;
	}
	if (StaticRestrictedMatrix) {
		for (int i = 0; i < DimOfProb; i++) {
			for (int j = 0; j < DimOfProb; j++) {
				clear_dof_matrix(StaticRestrictedMatrix[i][j]);
			}
			free(StaticRestrictedMatrix[i]);
		}
		free(StaticRestrictedMatrix);
		StaticRestrictedMatrix = NULL;
	}
	if (Fh) {
		for (int i = 0; i < DimOfProb; i++)
			free_dof_real_vec(Fh[i]);
		free(Fh);
		Fh = NULL;
	}
	if (Uh) {
		for (int i = 0; i < DimOfProb; i++)
			free_dof_real_vec(Uh[i]);
		free(Uh);
		Uh = NULL;
	}
	if (ProblemMatrix) {
		free_dof_matrix(ProblemMatrix);
		ProblemMatrix=NULL;
	}
	if (ProblemRHS) {
		free_dof_real_vec(ProblemRHS);
		ProblemRHS = NULL;
	}
	if (ProblemSol) {
		free_dof_real_vec(ProblemSol);
		ProblemSol = NULL;
	}
	if (pRestrictedVertexIndex) {
		for (int i = 0; i < DimOfProb; i++)
			free_dof_real_vec(pRestrictedVertexIndex[i]);
		free(pRestrictedVertexIndex);
		pRestrictedVertexIndex = NULL;
	}
	if (pExtendedVertexIndex) {
		for (int i = 0; i < DimOfProb; i++)
			free_dof_real_vec(pExtendedVertexIndex[i]);
		free(pExtendedVertexIndex);
		pExtendedVertexIndex = NULL;
	}
	if (IdxOfDofs) {
		delete []IdxOfDofs;
		IdxOfDofs=NULL;
	}
	if (NumOfBasFct) {
		free(NumOfBasFct);
		NumOfBasFct = NULL;
	}

	if (pRestrictedMesh) {
		free_mesh(pRestrictedMesh);
		pRestrictedMesh = NULL;
	}
	if (RestrictedFEspace) { /* not to be set before pRestrictedMesh */
		/* cleared by free_mesh(pRestrictedMesh)
		 for (int i = 0; i < DimOfProb; i++)
		 free_fe_space((FE_SPACE*)RestrictedFEspace[i]);
		 */
		free(RestrictedFEspace);
		RestrictedFEspace = NULL;
	}

	if (pMesh) {
		free_mesh(pMesh);
		pMesh = NULL;
	}
	if (FEspace) { /* not to be set before pMesh */
		/* cleared by free_mesh(pMesh)
		 for (int i = 0; i < DimOfProb; i++)
		 free_fe_space((FE_SPACE*)FEspace[i]);
		 */
		free(FEspace);
		FEspace = NULL;
	}
	if (ProblemFEspace) { /* not to be set before pMesh */
		/* cleared by free_mesh(pMesh) 
		 free_fe_space((FE_SPACE*)ProblemFEspace);
		 */
		ProblemFEspace = NULL;
	}
	if (pBoundaryEl) {
		for (int i = 0; i < DimOfWorld; i++) {
			free(pBoundaryEl[i][0]);
			free(pBoundaryEl[i][1]);
			free(pBoundaryEl[i]);
		}
		free(pBoundaryEl);
		pBoundaryEl = NULL;
	}
	if (pAllBoundaryEl) {
		free(pAllBoundaryEl);
		pAllBoundaryEl = NULL;
	}
	if (ConstrMatrix) {
		for (int i = 0; i < DimOfProb; i++) {
			free_dof_matrix(ConstrMatrix[i]);
		}
		free(ConstrMatrix);
		ConstrMatrix = NULL;
	}
	if (ConstrMatrixT) {
		for (int i = 0; i < DimOfProb; i++) {
			free_dof_matrix(ConstrMatrixT[i]);
		}
		free(ConstrMatrixT);
		ConstrMatrixT = NULL;
	}
	if (pOppositeDOFs) {
		free_dof_int_vec(pOppositeDOFs);
		pOppositeDOFs = NULL;
	}
}

/**
 Meshes pMesh and pRestrictedMesh are scaled by ScaleFactor
 and shifted by dx-ScaleFactor/2.0
 */
void CAlbertaInterface::TransformMeshes(REAL* dx/**center of transformed domain*/) {
	int i, j, k;
	for (int i = 0; i<DimOfWorld; i++)
		dx[i] -= ScaleFactor/2.0;

	for (i = 0; i<pMesh->n_elements; i++)
		for (j = 0; j<N_VERTICES(DimOfWorld); j++)
			for (k = 0; k<DimOfWorld; k++) {
				switch (k) {
				case 0:
					pMesh->macro_els[i].coord[j][k] = pVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k] *ScaleFactor
							+dx[0];
					break;
				case 1:
					pMesh->macro_els[i].coord[j][k] = pVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k] *ScaleFactor
							+dx[1];
					break;
				case 2:
					pMesh->macro_els[i].coord[j][k] = pVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k] *ScaleFactor
							+dx[2];
					break;
				}
			}

	for (i = 0; i<pRestrictedMesh->n_elements; i++)
		for (j = 0; j<N_VERTICES(DimOfWorld); j++)
			for (k = 0; k<DimOfWorld; k++) {
				switch (k) {
				case 0:
					pRestrictedMesh->macro_els[i].coord[j][k]
							= pRestrictedVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k]
									*ScaleFactor +dx[0];
					break;
				case 1:
					pRestrictedMesh->macro_els[i].coord[j][k]
							= pRestrictedVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k]
									*ScaleFactor +dx[1];
					break;
				case 2:
					pRestrictedMesh->macro_els[i].coord[j][k]
							= pRestrictedVertices[(i*N_VERTICES(DimOfWorld)+j)*DimOfWorld+k]
									*ScaleFactor +dx[2];
					break;
				}
			}
}

/**
 Local matrix cleanup
 */
void CAlbertaInterface::ClearLocalMatrix() {
	for (int i = 0; i < DimOfProb; i++)
		for (int j = 0; j < DimOfProb; j++)
			for (int k = 0; k < NumOfBasFct[i]; k++)
				for (int l = 0; l < NumOfBasFct[j]; l++)
					LocalMatrix[i][j][k][l] = 0;
}

/**
 Function fills pRestrictedVertexIndex what is array mapping elements
 from inside of domain (restricted domain) to the whole domain.
 */
int CAlbertaInterface::FillRestrictedVertexIndex() {
	TRAVERSE_STACK *Stack;
	TRAVERSE_STACK *StackRes;
	const DOF_ADMIN *Admin;
	const DOF_ADMIN *AdminRes;
	const EL_INFO *El_info;
	const EL_INFO *El_infoRes;
	DOF* pDofs= NULL;
	DOF* pDofsRes= NULL;

	for (int i = 0; i < DimOfProb; i++) {
		pDofs = (DOF*)realloc(pDofs, NumOfBasFct[i]*sizeof(DOF));
		pDofsRes = (DOF*)realloc(pDofsRes, NumOfBasFct[i]*sizeof(DOF));
		for (int j = 0; j < pExtendedVertexIndex[i]->size; j++)
			pExtendedVertexIndex[i]->vec[j] = -1.0;

		Admin = FEspace[i]->admin;
		AdminRes = RestrictedFEspace[i]->admin;
		dof_set(-1, pRestrictedVertexIndex[i]);
		Stack = get_traverse_stack();
		StackRes = get_traverse_stack();
		El_infoRes = traverse_first(StackRes, pRestrictedMesh, -1, 
		CALL_LEAF_EL|FILL_COORDS);
		for (El_info = traverse_first(Stack, pMesh, -1, 
		CALL_LEAF_EL|FILL_COORDS); El_info; El_info = traverse_next(Stack,
				El_info)) {
			if (!IsInside(El_info))
				continue;
			FEspace[i]->bas_fcts->get_dof_indices(El_info->el, Admin, pDofs);
			RestrictedFEspace[i]->bas_fcts->get_dof_indices(El_infoRes->el,
					AdminRes, pDofsRes);
			for (int j = 0; j < NumOfBasFct[i]; j++) {
				pRestrictedVertexIndex[i]->vec[pDofsRes[j]] = pDofs[j];
				pExtendedVertexIndex[i]->vec[pDofs[j]]=pDofsRes[j];
			}
			El_infoRes = traverse_next(StackRes, El_infoRes);
		}
	}

	if (pDofs) {
		free(pDofs);
		pDofs = NULL;
	}
	if (pDofsRes) {
		free(pDofsRes);
		pDofsRes = NULL;
	}
	return 0;
}

/**
 Function sets all dirichlet boundary in Matrix to the value Val.
 */
void CAlbertaInterface::PostProcessMatrix(REAL Val/** Value to be set on dirichlet boundary*/) {
	TRAVERSE_STACK *stack = get_traverse_stack();
	const EL_INFO *el_info;
	const S_CHAR *boundary;
	MATRIX_ROW *row;
	MATRIX_ROW *row_next;

	el_info = traverse_first(stack, pMesh, -1, CALL_LEAF_EL
	| FILL_BOUND | FILL_COORDS);
	dof_set(0.0, ProblemRHS);

	while (el_info) {
		for (int i = 0; i < DimOfProb; i++) {
			boundary = GetBound[i](el_info, nil);
			GetGlobDof[i](el_info->el, FEspace[i]->admin, Dofs[i]);
			for (int j = 0; j < NumOfBasFct[i]; j++) {
				if (boundary[j] == DIRICHLET) {
					row = ProblemMatrix->matrix_row[GetProblemDof(i,Dofs[i][j])];
					if (ProblemRHS->vec[GetProblemDof(i,Dofs[i][j])] == 0) {
						while (row) {
							row_next = row->next;
							free_matrix_row(ProblemMatrix->row_fe_space, row);
							row = row_next;
						}
						row = get_matrix_row(ProblemMatrix->row_fe_space);
						row->col[0] = GetProblemDof(i, Dofs[i][j]);
						row->entry[0] = Val;
						ProblemMatrix->matrix_row[GetProblemDof(i,Dofs[i][j])] = row;
						ProblemRHS->vec[GetProblemDof(i,Dofs[i][j])]=1;
					}
				}
			}
		}
		el_info = traverse_next(stack, el_info);
	}

	free_traverse_stack(stack);
	/*S_CHAR  vertex_bound[N_VERTICES(DimOfWorld)];
	 S_CHAR  edge_bound[N_EDGES_MAX];
	 S_CHAR  face_bound[N_FACES_MAX];*/
}

void CAlbertaInterface::AssembleProblemMatrix(const int **pMask,
		const DOF_MATRIX ***pMatrix) {
	MATRIX_ROW* row;
	MATRIX_ROW* rowExt;
	MATRIX_ROW* rowExtOld;
	int jcol;

	clear_dof_matrix(ProblemMatrix);

	for (int i = 0; i < DimOfProb; i++) {
		for (int j = 0; j < DimOfProb; j++) {
			if (pMask[i*DimOfProb+j])
				for (int dof = 0; dof < FEspace[i]->admin->size_used; dof++) {
					rowExtOld = ProblemMatrix->matrix_row[GetProblemDof(i, dof)];
					if (rowExtOld)
						while (rowExtOld->next)
							rowExtOld = rowExtOld->next;

					for (row = pMatrix[i][j]->matrix_row[dof]; row; row
							= row->next) {
						rowExt = get_matrix_row(ProblemMatrix->row_fe_space);
						for (int k = 0; k < ROW_LENGTH; k++) {
							jcol = row->col[k];
							if (ENTRY_USED(jcol)) {
								rowExt->col[k] = GetProblemDof(j, jcol);
								rowExt->entry[k] = row->entry[k];
							} else {
								rowExt->col[k] = UNUSED_ENTRY;
							}
						}
						if (!rowExtOld)
							ProblemMatrix->matrix_row[GetProblemDof(i, dof)] = rowExt;
						else
							rowExtOld->next = rowExt;
						rowExtOld = rowExt;
					}
				}
		}
	}
}

/** 
 Function extending weak formulation from elements inside domain into
 whole domain. Elements which are not inside are set to 0
 */
void CAlbertaInterface::ExtendRestrictedMatrix(int i, int j,
		const DOF_MATRIX *** pMatrix) {
	MATRIX_ROW* row;
	MATRIX_ROW* rowExt;
	MATRIX_ROW* rowExtOld;
	int jcol;

	clear_dof_matrix(ProblemMatrix);

	for (int dof = 0; dof < RestrictedFEspace[i]->admin->size_used; dof++) {
		rowExtOld = ProblemMatrix->matrix_row[GetProblemDof(i, (int)pRestrictedVertexIndex[i]->vec[dof])];
		if (rowExtOld)
			while (rowExtOld->next)
				rowExtOld = rowExtOld->next;
		for (row = pMatrix[i][j]->matrix_row[dof]; row; row = row->next) {
			rowExt = get_matrix_row(ProblemMatrix->row_fe_space);
			for (int k = 0; k < ROW_LENGTH; k++) {
				jcol = row->col[k];
				if (ENTRY_USED(jcol)) {
					rowExt->col[k] = GetProblemDof(j,
							(int)pRestrictedVertexIndex[i]->vec[jcol]);
					rowExt->entry[k] = row->entry[k];
				} else {
					rowExt->col[k] = UNUSED_ENTRY;
				}
			}
			if (!rowExtOld)
				ProblemMatrix->matrix_row[GetProblemDof(i, (int)pRestrictedVertexIndex[i]->vec[dof])]
						= rowExt;
			else
				rowExtOld->next = rowExt;
			rowExtOld = rowExt;

		}
	}
}

static inline REAL Length(const REAL_D x) {
	REAL sum = 0;
	for (int i = 0; i< DIM_OF_WORLD; i++) {
		sum += x[i]*x[i];
	}
	return sqrt(sum);
}

static inline REAL Abs(REAL x) {
	return x > 0 ? x : -x;
}

static inline REAL *CrossProd(REAL *v1, REAL *v2, REAL *out) {
	static REAL_D val;
	REAL *prd = out ? out : val;

	prd[0] = v1[1]*v2[2] - v1[2]*v2[1];
	prd[1] = v1[2]*v2[0] - v1[0]*v2[2];
	prd[2] = v1[0]*v2[1] - v1[1]*v2[0];

	return prd;
}

REAL CAlbertaInterface::GetInnerRadius1(const EL_INFO* el_info) {
	REAL line;
	line = el_info->coord[1][0] - el_info->coord[0][0];
	return Abs(line/2);
}

REAL CAlbertaInterface::GetInnerRadius2(const EL_INFO* el_info) {
	REAL_D line[3];
	REAL Sum = 0;
	for (int i = 0; i < DIM_OF_WORLD; i++) {
		line[0][i] = el_info->coord[1][i] - el_info->coord[0][i];
		line[1][i] = el_info->coord[2][i] - el_info->coord[0][i];
		line[2][i] = el_info->coord[2][i] - el_info->coord[1][i];
	};

	for (int i = 0; i < DIM_OF_WORLD; i++) {
		Sum += Length(line[i]);
	};

	return 2*el_volume(el_info)/(Sum*sqrt(2));
}

REAL CAlbertaInterface::GetInnerRadius3(const EL_INFO* el_info) {
	REAL_D line[4];
	REAL Surface[4];
	REAL Sur = 0;
	for (int i = 0; i < DIM_OF_WORLD; i++) {
		line[0][i] = el_info->coord[1][i]-el_info->coord[0][i]; //0
		line[1][i] = el_info->coord[3][i]-el_info->coord[0][i]; //2
		line[2][i] = el_info->coord[2][i]-el_info->coord[1][i]; //3
		line[3][i] = el_info->coord[3][i]-el_info->coord[2][i]; //5
	};

	Surface[0] = Length(CrossProd(line[0], line[1], NULL))/2.;
	Surface[1] = Length(CrossProd(line[1], line[3], NULL))/2.;
	Surface[2] = Length(CrossProd(line[3], line[2], NULL))/2.;
	Surface[3] = Length(CrossProd(line[2], line[0], NULL))/2.;

	for (int i = 0; i < DIM_OF_WORLD + 1; i++)
		Sur += Surface[i];

	return 6*el_volume(el_info)/(Sur*sqrt(3));

}

REAL CAlbertaInterface::GetInnerRadius(const EL_INFO* el_info) {
	REAL RetVal = 0;
	switch (DimOfWorld) {
	case 1:
		RetVal = GetInnerRadius1(el_info);
		break;
	case 2:
		RetVal = GetInnerRadius2(el_info);
		break;
	case 3:
		RetVal = GetInnerRadius3(el_info);
		break;
	}
	return RetVal;
}

REAL CAlbertaInterface::ComputeMeshVolume(MESH* mesh) {

	TRAVERSE_STACK *stack;
	const EL_INFO *el_info;
	REAL sum = 0;

	stack = get_traverse_stack();

	for (el_info = traverse_first(stack, mesh, -1, 
	CALL_LEAF_EL|FILL_COORDS); el_info; el_info = traverse_next(stack, el_info))
		sum += el_volume(el_info);

	return sum;
}

void CAlbertaInterface::SumMatrices(const int **pMask) {
	if (StaticMatrix)
		SumSystemMatrices(Matrix, StaticMatrix, pMask);
}

void CAlbertaInterface::SumRestrictedMatrices(const int **pMask) {
	if (StaticRestrictedMatrix)
		SumSystemMatrices(RestrictedMatrix, StaticRestrictedMatrix, pMask);
}

void CAlbertaInterface::SumSystemMatrices(DOF_MATRIX ***pDst,
		DOF_MATRIX ***pSrc, const int **pMask) {
	static REAL **pStatMat= NULL;
	if (!pStatMat) {
		pStatMat = (REAL**)malloc(sizeof(REAL*));
		pStatMat[0] = (REAL*)malloc(sizeof(REAL));
	}
	MATRIX_ROW* row;
	int jcol = 0;
	for (int i = 0; i<DimOfProb; i++)
		for (int j = 0; j<DimOfProb; j++) {
			if (pMask[i*DimOfProb+j]) {
				for (int dof = 0; dof
						< pSrc[i][j]->row_fe_space->admin->size_used; dof++) {
					for (row = pSrc[i][j]->matrix_row[dof]; row; row
							= row->next) {
						for (int k = 0; k<ROW_LENGTH; k++) {
							jcol = row->col[k];
							if (ENTRY_USED(jcol))
								pStatMat[0][0] = row->entry[k];
							add_element_matrix(pDst[i][j], 1., 1, 1, &dof,
									&jcol, (const REAL**)pStatMat, NULL);
						}
					}
				}
			}
		}
}

void CAlbertaInterface::ClearMatrix() {
	for (int i = 0; i < DimOfProb; i++)
		for (int j = 0; j < DimOfProb; j++)
			clear_dof_matrix(Matrix[i][j]);
}

void CAlbertaInterface::ClearRestrictedMatrix() {
	for (int i = 0; i < DimOfProb; i++)
		for (int j = 0; j < DimOfProb; j++)
			clear_dof_matrix(RestrictedMatrix[i][j]);
}

void CAlbertaInterface::GetNumOfBoundEl() {
	TRAVERSE_STACK *stack;
	const EL_INFO *el_info;
	int Counter = 0;
	int CounterIn = 0;

	stack = get_traverse_stack();

	for (el_info = traverse_first(stack, pMesh, -1, 
	CALL_LEAF_EL|FILL_BOUND|FILL_COORDS); el_info; el_info = traverse_next(
			stack, el_info)) {
		CounterIn=0;
		for (int j=0; j<DimOfWorld+1; j++)
			if (el_info->coord[j][0]==0)
				CounterIn++;
		if (CounterIn==DimOfWorld)
			Counter++;

	}
	free_traverse_stack(stack);
	NumOfBoundEl=Counter;
}

int CAlbertaInterface::AreCorrespondingEls(MACRO_EL* macro_el, REAL * y,
		int orientation) {
	if (DimOfWorld == 1)
		return 1;
	REAL x[2] = { 0 };
	ComputeCenterOfGravity(macro_el, x, orientation, 0);
	for (int i = 0; i < DimOfWorld-1; i++) {
		if (x[i]!=y[i])
			return 0;
	}
	return 1;
}

void CAlbertaInterface::CoincidentBoundToEl(const EL_INFO* el_info, REAL* x) {
	static REAL xx[3];
	x=x ? x : xx;
	for (int i=0; i<3; i++) {
		x[i]=-1;
	}

	int Counter=0;
	for (int i=0; i<DimOfWorld; i++) {
		Counter=0;
		for (int j=0; j<DimOfWorld+1; j++)
			if (el_info->coord[j][i]==0)
				Counter++;
		if (Counter==DimOfWorld)
			x[i]=1;
	}
}

void CAlbertaInterface::CoincidentOppositeBoundToEl(const EL_INFO* el_info,
		REAL* x) {

	for (int i=0; i<3; i++) {
		x[i]=-1;
	}

	int Counter=0;
	for (int i=0; i<DimOfWorld; i++) {
		Counter=0;
		for (int j=0; j<DimOfWorld+1; j++)
			if (el_info->coord[j][i]==1)
				Counter++;
		if (Counter==DimOfWorld)
			x[i] = 1;
	}
}

void CAlbertaInterface::FillPrimaryBoundElArray() {
	TRAVERSE_STACK *stack;
	const EL_INFO *el_info;
	REAL x[3]= { 0 };
	int Counter[3] = { 0 };

	stack = get_traverse_stack();

	for (el_info = traverse_first(stack, pMesh, -1, 
	CALL_LEAF_EL|FILL_BOUND|FILL_COORDS); el_info; el_info = traverse_next(
			stack, el_info)) {
		CoincidentBoundToEl(el_info, x);
		for (int i=0; i<DimOfWorld; i++) {
			if (x[i]!=-1)
				pBoundaryEl[i][0][Counter[i]++]=(MACRO_EL*)el_info->macro_el;
		}
	}
	free_traverse_stack(stack);
}

void CAlbertaInterface::ComputeCenterOfGravity(MACRO_EL *macro_el,
		REAL * center, int orientation, int par) {
	if (DimOfWorld == 1)
		return;
	int Index;
	for (int i = 0; i < DimOfWorld+1; i++)
		if (macro_el->coord[i][orientation] != par) {
			Index=0;
			for (int j = 0; j < DimOfWorld; j++) {
				if (j!=orientation)
					center[Index++]+=macro_el->coord[i][j];
			}
		}
}

void CAlbertaInterface::FillOppositeDOFs(MACRO_EL *pFrontEl, MACRO_EL *pRearEl,
		int orientation) {
	static const DOF *(*get_dof)(const EL *, const DOF_ADMIN *, DOF *);
	const DOF *pRearDof= NULL;
	EL_INFO FrontInfo, RearInfo;
	REAL_D FrontPoint, RearPoint;
	double *pFront, *pRear;
	int i, j, k, opposite;

	/* create minimal EL_INFO */
	FrontInfo.mesh = pMesh;
	FrontInfo.fill_flag = FILL_COORDS;
	RearInfo.mesh = pMesh;
	RearInfo.fill_flag = FILL_COORDS;
	for (int i = 0; i < N_VERTICES(DimOfWorld); i++) {
		for (int j = 0; j < DimOfWorld; j++) {
			(FrontInfo.coord[i])[j]=(pFrontEl->coord[i])[j];
			(RearInfo.coord[i])[j]=(pRearEl->coord[i])[j];
		}
	}

	get_dof = LinearFEspace->bas_fcts->get_dof_indices;

	/* get dofs on front and rear face */
	pFrontDof = get_dof(pFrontEl->el, LinearFEspace->admin, (DOF*)pFrontDof);
	pRearDof = get_dof(pRearEl->el, LinearFEspace->admin, nil);

	for (i = 0; i < LinearFEspace->bas_fcts->n_bas_fcts; i++) {
		pFront=(double*)coord_to_world((const EL_INFO*)&FrontInfo, (LAGRANGE_NODES(LinearFEspace->bas_fcts))[i], FrontPoint);
		if (pFront[orientation] != 0.) // not a boundary element
			continue;
		for (j = 0; j < LinearFEspace->bas_fcts->n_bas_fcts; j++) {
			pRear=(double*)coord_to_world((const EL_INFO*)&RearInfo, (LAGRANGE_NODES(LinearFEspace->bas_fcts))[j], RearPoint);
			if (pRear[orientation] != 1.) // not a boundary element
				continue;
			opposite = 1;
			for (k = 0; k < DimOfWorld; k++) {
				if ((pFront[k] != pRear[k]) && (k!=orientation)) {
					opposite = 0;
					break;
				}
			}
			if (opposite) {
				pOppositeDOFs->vec[pFrontDof[i]]=pRearDof[j];
				pOppositeDOFs->vec[pRearDof[j]]=pFrontDof[i];
				break;
			}
		}
	}
}

void CAlbertaInterface::FillOppositeBoundElArray() {
	TRAVERSE_STACK *stack;
	const EL_INFO *el_info;
	REAL x[3] = { 0 };
	REAL y[2] = { 0 };

	stack = get_traverse_stack();

	pFrontDof = new DOF[LinearFEspace->bas_fcts->n_bas_fcts]; // temporary for fct FillOppositeDOFs

	for (el_info = traverse_first(stack, pMesh, -1, 
	CALL_LEAF_EL|FILL_BOUND|FILL_COORDS); el_info; el_info = traverse_next(
			stack, el_info)) {
		CoincidentOppositeBoundToEl(el_info, x);
		for (int i = 0; i<DimOfWorld; i++) {
			if (x[i] != -1) {
				y[0] = y[1] = 0.0;
				ComputeCenterOfGravity((MACRO_EL*)el_info->macro_el, y, i, 1);
				for (int j = 0; j < NumOfBoundEl; j++) {
					if (AreCorrespondingEls(pBoundaryEl[i][0][j], y, i)) {
						FillOppositeDOFs(pBoundaryEl[i][0][j],
								(MACRO_EL*)el_info->macro_el, i);
						pBoundaryEl[i][1][j] = (MACRO_EL*)el_info->macro_el;
						break;
					}
				}
			}
		}
	}
	free_traverse_stack(stack);

	if (pFrontDof) {
		delete pFrontDof;
		pFrontDof = NULL;
	}
}

void CAlbertaInterface::AssembleConstrMatrices() {
#if 0
	for (int i = 0; i< DimOfProb; i++) {
		matrix_info = nil;
		OPERATOR_INFO o_info = {nil};
		o_info.row_fe_space = LinearFEspace;
		o_info.col_fe_space = FEspace[i];
		o_info.init_element = init_element;
		o_info.c = c;
		//o_info.quad           = quad;
		//o_info.use_get_bound = true; // Dirichlet boundary conditions! 
		//o_info.user_data = MEM_ALLOC(1, struct op_info); // user data!
		o_info.fill_flag = CALL_LEAF_EL|FILL_COORDS;

		matrix_info = fill_matrix_info(&o_info, nil);
		fill = matrix_info->el_matrix_fct;
		loc_mat = fill(el_info, matrix_info->fill_info);

		clear_dof_matrix(ConstrMatrix[i]); // assembling of matrix
		update_matrix(ConstrMatrix[i], matrix_info);

	}

	for (int i = 0; i< DimOfProb; i++) {
		for (int j = 0; j< NumOfBoundEl; j++) {

		}
	}
#endif
}

