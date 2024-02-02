#include "SaveToVTK.h"

// define taken from macro.c
#define VERT_IND(dim,i,j) ((i)*N_VERTICES(dim)+(j)) 

using namespace std;

int SaveToVTK::ActComponentIdx = 0;

void SaveToVTK::FillDofIndices(DOF_INT_VEC * pDst, MESH * pMesh) {
	TRAVERSE_STACK *stack;
	const DOF_ADMIN *admin;
	//	FE_SPACE fe_space = { "write fe_space", nil, nil };
	const EL_INFO *el_info;
	DOF_INT_VEC *dof_vert_ind;
	int dim = pMesh->dim, n0, nv, i, *vert_ind= nil;

	admin = get_vertex_admin(pMesh);

	n0 = admin->n0_dof[VERTEX];
	//	fe_space.admin = admin;

	dof_vert_ind = get_dof_int_vec("vertex indices", pDst->fe_space);
	GET_DOF_VEC(vert_ind, dof_vert_ind);
	FOR_ALL_DOFS(admin, vert_ind[dof] = -1);

	stack = get_traverse_stack();

	nv = 0;
	FOR_ALL_DOFS(admin, vert_ind[dof] = -1);

	/* The first pass assigns global indices to the vertices */
	for (el_info = traverse_first(stack, pMesh, -1, CALL_LEAF_EL | FILL_COORDS); el_info; el_info
			= traverse_next(stack, el_info)) {
		for (i = 0; i < N_VERTICES(dim); i++)
			if (vert_ind[el_info->el->dof[i][n0]] == -1) {
				/* assign a global index to each vertex                                     */
				vert_ind[el_info->el->dof[i][n0]] = nv;
				nv++;
			}
	}
	/* The second pass assigns vert_indices to pDst */
	// maybe this can be done also in the previous cycle...
	for (el_info = traverse_first(stack, pMesh, -1, CALL_LEAF_EL | FILL_COORDS); el_info; el_info
			= traverse_next(stack, el_info)) {
		for (i = 0; i < N_VERTICES(dim); i++)
			pDst->vec[vert_ind[el_info->el->dof[i][n0]]]
					= el_info->el->dof[i][n0];
	}

	free_dof_int_vec(dof_vert_ind);
	free_traverse_stack(stack);
}

SaveToVTK::SaveToVTK(EXACT_SOL exact, int world, int prob) {
	pDofInd = NULL;
	pMesh = NULL;
	pMacroData = NULL;
	DimOfProb = prob;
	DimOfWorld = world;
	if (exact)
		pExact = exact;
}

REAL SaveToVTK::LocalExact(const REAL_D x) {
	return (pExact) ? pExact(x, ActComponentIdx) : 0;
}

void SaveToVTK::SaveVTKfile(const DOF_REAL_VEC **pSol, const char *pName) {
	FUNCNAME("SaveToVTK::SaveVTKfile");
	
	FILE * fpOut;
	REAL_D x;
	int i, j;

	if (!pSol) {
		ERROR("Solution missing. Possible problem in the library.\n");
		exit(ERC_SOLUTION_MISSING);
	}
	if (!pName){
		WARNING("File name not set. Cannot continue in file save.\n");
		return;
	}

	if (pMesh != pSol[0]->fe_space->mesh) { // mesh has changed
		if (pMacroData)
			free_macro_data(pMacroData);
		if (pDofInd)
			free_dof_int_vec(pDofInd);
		pMesh=pSol[0]->fe_space->mesh;
		pDofInd = get_dof_int_vec("pDofInd", pSol[0]->fe_space);
		if (!(pMacroData = mesh2macro_data(pMesh))) {
			printf("mesh2macro error\n");
			return;
		}
		FillDofIndices(pDofInd, pMesh);
	}

	fpOut = fopen(pName, "w");

	fprintf(
			fpOut,
			"#0 This is not a valid VTK file. In order to create such, all lines contaning #X\n");
	fprintf(
			fpOut,
			"#0 where X is a integer should be consulted and removed or changed in desired\n");
	fprintf(fpOut, "#0 way. Possible 'X' values:\n");
	fprintf(fpOut, "#0 0 - comment\n");
	fprintf(
			fpOut,
			"#0 1 - line containing problem description: DimOfProb XX DimOfWorld XX ExactSol XX\n");
	fprintf(fpOut,
			"#0     ExactSol is 1 if in results exact solution is provided, 0 otherwise\n");
	fprintf(
			fpOut,
			"#0 2 - Leave this lines (after removing starting sign #2) if output file should\n");
	fprintf(fpOut, "#0     contain scalar values.\n");
	fprintf(
			fpOut,
			"#0 3 - Leave this lines (after removing starting sign #3) if output file should\n");
	fprintf(fpOut,
			"#0     contain vector values. Vector has to have 3 components.\n");
	fprintf(
			fpOut,
			"#0 4 - DimOfProb floating point values of solution followed (if ExactSol = 1) by\n");
	fprintf(
			fpOut,
			"#0     DimOfProb floating point values of exact solution. Every row represents values\n");
	fprintf(fpOut, "#0     for one node.\n");
	fprintf(fpOut, "#1 DimOfProb %d DimOfWorld %d ExactSol %d\n", DimOfProb,
			DimOfWorld, (pExact) ? 1 : 0);
	fprintf(fpOut, "# vtk DataFile Version 3.0\n");
	fprintf(fpOut, "Voids\n");
	fprintf(fpOut, "ASCII\n");
	fprintf(fpOut, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fpOut, "POINTS %i float\n", pMacroData->n_total_vertices);

	// vertex coordinates
	for (i = 0; i < pMacroData->n_total_vertices; i++) {
#if DIM_OF_WORLD == 1
		fprintf(fpOut, "%lf %lf %lf\n", pMacroData->coords[i][0], 0., 0.);
#elif DIM_OF_WORLD == 2
		fprintf(fpOut, "%lf %lf %lf\n", pMacroData->coords[i][0],
				pMacroData->coords[i][1], 0.);
#elif DIM_OF_WORLD == 3
		fprintf(fpOut, "%lf %lf %lf\n", pMacroData->coords[i][0],
				pMacroData->coords[i][1], pMacroData->coords[i][2]);

#endif
	}

	// cell definition
#if DIM_OF_WORLD == 1
	fprintf(fpOut, "\nCELLS %d %d\n", pMacroData->n_macro_elements, 3
			*pMacroData->n_macro_elements);
#elif DIM_OF_WORLD == 2
	fprintf(fpOut, "\nCELLS %d %d\n", pMacroData->n_macro_elements, 4
			*pMacroData->n_macro_elements);
#elif DIM_OF_WORLD == 3
	fprintf(fpOut, "\nCELLS %d %d\n", pMacroData->n_macro_elements, 5
			*pMacroData->n_macro_elements);
#endif

	// cell elements
	for (i = 0; i < pMacroData->n_macro_elements; i++) {
#if DIM_OF_WORLD == 1
		fprintf(fpOut, "%d %d %d\n", 2, pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,0)], pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,1)]);
#elif DIM_OF_WORLD == 2
		fprintf(fpOut, "%d %d %d %d\n", 3, pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,0)], pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,1)], pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,2)]);
#elif DIM_OF_WORLD == 3
		fprintf(fpOut, "%d %d %d %d %d\n", 4,
				pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,0)],
				pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,1)],
				pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,2)],
				pMacroData->mel_vertices[VERT_IND(DimOfWorld,i,3)]);
#endif		
	}

	// cell types
	fprintf(fpOut, "\nCELL_TYPES %d\n", pMacroData->n_macro_elements);
	for (i = 0; i < pMacroData->n_macro_elements; i++) {
#if DIM_OF_WORLD == 1
		fprintf(fpOut, "%d ", 3); // line
#elif DIM_OF_WORLD == 2
		fprintf(fpOut, "%d ", 5); // triangle
#elif DIM_OF_WORLD == 1
		fprintf(fpOut, "%d ", 10); // tetrahedron
#endif
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\nPOINT_DATA %d\n", pMacroData->n_total_vertices);
	fprintf(fpOut, "#2 SCALARS u_h float 1\n");
	fprintf(fpOut, "#2 LOOKUP_TABLE default\n");
	fprintf(fpOut, "#3 VECTORS m_h float\n");

	for (i = 0; i < pMacroData->n_total_vertices; i++) {
		fprintf(fpOut, "#4 ");
		for (j = 0; j < DimOfProb; j++)
			fprintf(fpOut, "%lf ", pSol[j]->vec[pDofInd->vec[i]]);
		if (pExact) {
			for (j = 0; j < DIM_OF_WORLD; j++)
				x[j] = pMacroData->coords[i][j];
			for (j = 0; j < DimOfProb; j++) {
				ActComponentIdx = j;
				fprintf(fpOut, "%lf ", LocalExact(x));
			}
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);
}

SaveToVTK::~SaveToVTK() {
	if (pMacroData) {
		free_macro_data(pMacroData);
		pMacroData = NULL;
	}
	if (pDofInd) {
		free_dof_int_vec(pDofInd);
		pDofInd = NULL;
	}
}

