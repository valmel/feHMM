#include <stdio.h>
#include <iostream>
#include <fstream>
#include "feHMM.h"
#include <math.h>
#include <gsl/gsl_integration.h>

CUserData MyData;

gsl_integration_workspace * w;
gsl_function F;
const double HomEps=0.05;
const double Pi=3.141592653589793238;

static int DimOfProb = 0;

#ifdef USING_MPI
int Rank;
#endif

const QUAD **quad= NULL;
const QUAD_FAST **quad_fast= NULL;
REAL ***phi= NULL;
int *n_phi= NULL;
REAL **wdetf_qp= NULL;

inline double A(double x, double y) {
	return cos(2*Pi*x)+cos(2*Pi*y)+3;
}

inline double Term1(double x, void * params) {
	return x/A(x/HomEps, x/HomEps*1./HomEps);
}

inline double Term2(double x, void * params) {
	return 1./A(x/HomEps, x/HomEps*1./HomEps);
}

using namespace std;

/*--------------------------------------------------------------------------*/
/* struct ellipt_leaf_data: structure for storing one REAL value on each    */
/*                          leaf element as LEAF_DATA                       */
/* rw_el_est():  return a pointer to the memory for storing the element     */
/*               estimate (stored as LEAF_DATA), called by ellipt_est()     */
/* get_el_est(): return the value of the element estimates (from LEAF_DATA),*/
/*               called by adapt_method_stat() and graphics()               */
/*--------------------------------------------------------------------------*/

struct ellipt_leaf_data {
	REAL estimate; /*  one real for the estimate                   */
};

/*--------------------------------------------------------------------------*/
/* For test purposes: exact solution and its gradient (optional)            */
/*--------------------------------------------------------------------------*/

static REAL u(const REAL_D x, int i) {
	static REAL RetVal;
	RetVal = 0;
	switch (MyData.ProblemNum) {
	case 1:
		RetVal = 1;
		break;
	case 2:
		RetVal = (exp(-10.0*SCP_DOW(x, x)));
		break;
	case 3:
		double a, b, c, d, error;
		F.function = &Term1;
		F.params = NULL;
		gsl_integration_qags(&F, 0, x[0], 0, 1e-8, 1000000, w, &a, &error);
		gsl_integration_qags(&F, 0, 1., 0, 1e-8, 1000000, w, &b, &error);
		F.function = &Term2;
		gsl_integration_qags(&F, 0, 1., 0, 1e-8, 1000000, w, &c, &error);
		gsl_integration_qags(&F, 0, x[0], 0, 1e-8, 1000000, w, &d, &error);
		RetVal = -a+b/c*d;
		break;
	}
	return RetVal;
}

/*static const REAL *grd_u(const REAL_D x, REAL_D input) {
 static REAL_D buffer = {};
 REAL *grd = input ? input : buffer;

 REAL          ux = exp(-10.0*SCP_DOW(x,x));
 int           n;

 for (n = 0;  n < DIM_OF_WORLD; n++)
 grd[n] = -20.0*x[n]*ux;

 return(grd);
 }*/

/*--------------------------------------------------------------------------*/
/* problem data: right hand side, boundary values                           */
/*--------------------------------------------------------------------------*/

static REAL g(const REAL_D x, int i) /* boundary values, not optional */
{
	switch (MyData.ProblemNum) {
	case 1:
		return (u(x, i));
	case 2:
		return (u(x, i));
	case 3:
		return 0;
	}
	return 0;
}

static REAL f(const REAL_D x, int i) /* -Delta u, not optional        */
{
	REAL r2, ux;
	switch (MyData.ProblemNum) {
	case 1:
		return 0.;
	case 2: {
		r2 = SCP_DOW(x, x);
		ux = exp(-10.0*r2);
		return (-(400.0*r2 - 20.0*DIM_OF_WORLD)*ux);
	}
	case 3:
		return 1.;
	}
	return 0;
}

struct op_info {
	REAL_D Lambda[N_LAMBDA]; /*  the gradient of the barycentric coordinates */
	REAL det; /*  |det D F_S|                                 */
};

static int init_element(const EL_INFO *el_info, const QUAD *quad[3], void *ud) {
	FUNCNAME("init_element");
	struct op_info *info = (struct op_info *)ud;

	switch (el_info->mesh->dim) {
	case 1:
		info->det = el_grd_lambda_1d(el_info, info->Lambda);
		break;
#if DIM_OF_WORLD > 1
	case 2:
		info->det = el_grd_lambda_2d(el_info, info->Lambda);
		break;
#if DIM_OF_WORLD > 2
		case 3:
		info->det = el_grd_lambda_3d(el_info, info->Lambda);
		break;
#endif
#endif
	default:
		ERROR_EXIT("Illegal dim!\n");
	}
	return 0;
}

const REAL (*LALt(const EL_INFO *el_info, const QUAD *quad,
				int iq, void *ud))[N_LAMBDA] {
	struct op_info *info = (struct op_info *)ud;
	int i, j, k, dim = el_info->mesh->dim;
	static REAL LALt[N_LAMBDA][N_LAMBDA];
	static REAL Bary[DIM_OF_WORLD+1]; // = {1./N_LAMBDA};
	for (i = 0; i < DIM_OF_WORLD+1; i++)
		Bary[i] = 1./(DIM_OF_WORLD+1);
	//static const REAL* pBary=Bary;
	static REAL x[DIM_OF_WORLD];

	coord_to_world(el_info, Bary, x);
	for (i = 0; i <= dim; i++)
	for (j = i; j <= dim; j++) {
		for (LALt[i][j] = k = 0; k < DIM_OF_WORLD; k++)
		LALt[i][j] += info->Lambda[i][k]*info->Lambda[j][k];
		switch (MyData.ProblemNum) {
			case 1:
			LALt[i][j] *= info->det; break;
			case 2:
			LALt[i][j] *= info->det; break;
			case 3:
			LALt[i][j] *= info->det*A(x[0]/HomEps, x[0]/(HomEps*HomEps)); break;
		}
		LALt[j][i] = LALt[i][j];
	}
	return((const REAL (*)[N_LAMBDA]) LALt);
}

void build(const FE_SPACE **fe_space, DOF_MATRIX*** matrix) {
	MESH* pMesh=fe_space[0]->mesh;
	dof_compress(pMesh);
	MSG("%d DOFs for %s\n", fe_space[0]->admin->size_used, fe_space[0]->name);

	static const EL_MATRIX_INFO *matrix_info= nil;
	
	for (int i = 0; i < DimOfProb; i++) {
		matrix_info = nil;
		
		OPERATOR_INFO o_info = { nil };
		
		o_info.row_fe_space = fe_space[i];
		o_info.col_fe_space = fe_space[i];
		o_info.init_element = init_element;
		o_info.LALt = LALt;
		//o_info.quad           = quad;
		o_info.LALt_pw_const = true; // pw const. assemblage is faster 
		o_info.LALt_symmetric = true; // symmetric assemblage is faster 
		o_info.use_get_bound = true; // Dirichlet boundary conditions! 
		o_info.user_data = MEM_ALLOC(1, struct op_info); // user data! 
		o_info.fill_flag = CALL_LEAF_EL|FILL_COORDS;

		matrix_info = fill_matrix_info(&o_info, nil);

		clear_dof_matrix(matrix[i][i]); // assembling of matrix
		update_matrix(matrix[i][i], matrix_info);
	}
}


void buildRHS(const EL_INFO *El_info, const FE_SPACE **fe_space, REAL **rhs) {
	REAL det;
	REAL_D x;
	MESH* pMesh = fe_space[0]->mesh;
	double val;
	int i, j, iq;

	if (!quad) {
		n_phi = (int*)malloc(DimOfProb*sizeof(int));
		quad = (const QUAD**)malloc(DimOfProb*sizeof(QUAD*));
		quad_fast = (const QUAD_FAST**)malloc(DimOfProb*sizeof(QUAD_FAST*));
		phi = (REAL***)malloc(DimOfProb*sizeof(REAL**));
		wdetf_qp = (REAL**)malloc(DimOfProb*sizeof(REAL*));
		for (i = 0; i < DimOfProb; i++) {
			quad[i] = get_quadrature(pMesh->dim, 2
					*fe_space[i]->bas_fcts->degree - 2);
			quad_fast[i] = get_quad_fast(fe_space[i]->bas_fcts, quad[i], 
			INIT_PHI);
			phi[i] = quad_fast[i]->phi;
			n_phi[i] = fe_space[i]->bas_fcts->n_bas_fcts;
			wdetf_qp[i] = (REAL*)malloc(quad[i]->n_points*sizeof(REAL));
		}
	}

	//L2scp_fct_bas(f, quad, rhs[0]);....

	det = el_det(El_info);

	for (i = 0; i < DimOfProb; i++) {
		for (iq = 0; iq < quad[i]->n_points; iq++) {
			coord_to_world(El_info, quad[i]->lambda[iq], x);
			wdetf_qp[i][iq] = quad[i]->w[iq]*det*(*f)(x, i);
		}

		for (j = 0; j < n_phi[i]; j++) {
			for (val = iq = 0; iq < quad[i]->n_points; iq++)
				val += phi[i][iq][j]*wdetf_qp[i][iq];
			rhs[i][j] += val;
		}
	}
}

/*--------------------------------------------------------------------------*/
/* main program                                                             */
/*--------------------------------------------------------------------------*/

int main(int argc, char **argv) {
#ifdef USING_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
#endif
	change_msg_out(fopen("msg.info", "wt")); //redirection of output to a extern file
	msg_info=10;
	static int pMaskOfZorro[2][2] = {{1,0},{0,1}};

	w = gsl_integration_workspace_alloc(1000000);

	MyData.LoadUserDataFromFile("INIT/init.dat");

	MyData.pBuildDynamic = &build; // user has to build the global matrix (it is used at the lowest level)
	MyData.pSetUpRHS = &buildRHS; // user has to build the macro RHS, i.e., the RHS of the highest level
	MyData.pExactSolution = &u;
	MyData.pDirBoundCond = &g;
	MyData.pMatrixMask = (int**)&(pMaskOfZorro[0][0]);

    DimOfProb = MyData.DimOfProb;

	feHMM* pMyfeHMM;// initialize feHMM with user data
	ofstream OUT;
	OUT.open("out.txt");
	OUT<<"# X Y Z"<<endl;
	//for(int i=4; i<5; i++){
		//for(int j=3; j<5; j++){
		//MyData.pNumOfRefin[0] = i*2;
	    //MyData.pNumOfRefin[1] =	j*2;
		pMyfeHMM = new feHMM(&MyData);
		pMyfeHMM->Solve(false);
		OUT<<"E="<<pMyfeHMM->L2Error()<<endl;
		pMyfeHMM->SaveSolutionToVTK("result.vt_");
		//OUT<<i*2<<" "<<j*2<<" "<<pMyfeHMM->L2Error()<<endl;
		delete pMyfeHMM;
		//}
		//OUT<<endl;
	//}
	
#ifdef USING_MPI
	if (Rank == 0)
#endif
	{
		//pMyfeHMM.SaveSolutionToVTK("result.vtk");
		cout<<"Skoncil som"<<endl;
	}
	gsl_integration_workspace_free(w);
	return 0;
}

