#include <stdio.h>
#include <iostream>
#include "feHMMpar.h"
#include <math.h>
#include <gsl/gsl_integration.h>

CUserDataPar MyData;

gsl_integration_workspace * w;
gsl_function F;
const double HomEps=0.01;
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

inline
double A(double x) {
	return cos(2*Pi*x)+2;
}

inline
double Term1(double x, void * params) {
	return x/A(x/HomEps);
}

inline
double Term2(double x, void * params) {
	return 1./A(x/HomEps);
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

static REAL ustatic(const REAL_D x, const REAL t, int i) {
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
		gsl_integration_qags(&F, 0, x[0], 0, 1e-5, 100000, w, &a, &error);
		gsl_integration_qags(&F, 0, 1., 0, 1e-5, 100000, w, &b, &error);
		F.function = &Term2;
		gsl_integration_qags(&F, 0, 1., 0, 1e-5, 100000, w, &c, &error);
		gsl_integration_qags(&F, 0, x[0], 0, 1e-5, 100000, w, &d, &error);
		RetVal = -a+b/c*d;
		break;
	}
	return RetVal;
}

static REAL u(const REAL_D x, const REAL t, int i) {
	static REAL RetVal;
	RetVal = 0;
	switch (MyData.ProblemNum) {
	case 1:
		RetVal = t * ustatic(x, t, i);
		break;
	case 2:
		RetVal = ustatic(x, t, i);
		break;
	case 3:
		RetVal = t * ustatic(x, t, i);
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

static double ActualTime;

static REAL g(const REAL_D x, REAL t, int i) /* boundary values, not optional */
{
	switch (MyData.ProblemNum) {
	case 1:
		return (u(x, t, i));
	case 2:
		return (u(x, t, i));
	case 3:
		return 0;
	}
	return 0;
}

static REAL f(const REAL_D x, const REAL t, int i) /* -Delta u, not optional        */
{
	REAL r2, ux;
	switch (MyData.ProblemNum) {
	case 1:
		return ustatic(x, t, i);
	case 2: {
		r2 = SCP_DOW(x, x);
		ux = exp(-10.0*r2);
		return (-(400.0*r2 - 20.0*DIM_OF_WORLD)*ux);
	}
	case 3:
		return ustatic(x, t, i) + t;
	}
	return 0;
}

/*--------------------------------------------------------------------------*/
/* build(): assemblage of the linear system: matrix, load vector,           */
/*          boundary values, called by adapt_method_stat()                  */
/*          on the first call initialize u_h, f_h, matrix and information   */
/*          for assembling the system matrix                                */
/*                                                                          */
/* struct op_info: structure for passing information from init_element() to */
/*                 LALt()                                                   */
/* init_element(): initialization on the element; calculates the            */
/*                 coordinates and |det DF_S| used by LALt; passes these    */
/*                 values to LALt via user_data,                            */
/*                 called on each element by update_matrix()                */
/* LALt():         implementation of -Lambda id Lambda^t for -Delta u,      */
/*                 called by update_matrix() after init_element()           */
/*--------------------------------------------------------------------------*/

struct op_info {
	REAL_D Lambda[N_LAMBDA]; /*  the gradient of the barycentric coordinates */
	REAL det; /*  |det D F_S|                                 */
	REAL tau;
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
			LALt[i][j] *= info->det*A(x[0]/HomEps); break;
		}
		LALt[j][i] = LALt[i][j];
	}
	return((const REAL (*)[N_LAMBDA]) LALt);
}

static REAL c(const EL_INFO *el_info, const QUAD *quad, int iq, void *ud) {
	struct op_info *info = (struct op_info *)ud;
	//  DEBUG_TEST_EXIT(info->quad_fast->quad == quad, "quads differ\n");
	return (info->det/info->tau);
}

void build(const FE_SPACE **fe_space, DOF_MATRIX*** matrix, OLD_SOL_PAR uhold,
		REAL t, REAL tau) {
	MESH* pMesh=fe_space[0]->mesh;
	dof_compress(pMesh);
	MSG("%d DOFs for %s\n", fe_space[0]->admin->size_used, fe_space[0]->name);

	static const EL_MATRIX_INFO *matrix_info1= nil;
	static const EL_MATRIX_INFO *matrix_info2= nil;
	//if (matrix_info[i] && matrix_info[i]->row_admin->mesh!=pMesh)	

	for (int i = 0; i < DimOfProb; i++) {
		matrix_info1 = nil;
		matrix_info2 = nil;

		OPERATOR_INFO o_info1 = { nil };
		OPERATOR_INFO o_info2 = { nil };

		o_info1.row_fe_space = fe_space[i];
		o_info1.col_fe_space = fe_space[i];
		o_info1.init_element = init_element;
		o_info1.LALt = LALt;
		//o_info.quad           = quad;
		o_info1.LALt_pw_const = true; // pw const. assemblage is faster 
		o_info1.LALt_symmetric = true; // symmetric assemblage is faster 
		o_info1.use_get_bound = true; // Dirichlet boundary conditions! 
		o_info1.user_data = MEM_ALLOC(1, struct op_info); // user data! 
		o_info1.fill_flag = CALL_LEAF_EL|FILL_COORDS;

		matrix_info1 = fill_matrix_info(&o_info1, nil);

		o_info2.row_fe_space = fe_space[i];
		o_info2.col_fe_space = fe_space[i];
		o_info2.init_element = init_element;
		o_info2.c = c;
		//o_info.quad           = quad;
		o_info2.use_get_bound = true; // Dirichlet boundary conditions! 
		o_info2.user_data = MEM_ALLOC(1, struct op_info); // user data!
		((op_info*)o_info2.user_data)->tau = tau;
		o_info2.fill_flag = CALL_LEAF_EL|FILL_COORDS;

		matrix_info2 = fill_matrix_info(&o_info2, nil);

		clear_dof_matrix(matrix[i][i]); // assembling of matrix
		update_matrix(matrix[i][i], matrix_info1);
		update_matrix(matrix[i][i], matrix_info2);
	}
}

void buildRHS(const EL_INFO *El_info, const FE_SPACE **fe_space, REAL **rhs,
		OLD_SOL_PAR uhold, REAL t, REAL tau) {

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
			wdetf_qp[i][iq] = quad[i]->w[iq]*det*((*f)(x, t, i)+1./tau*(uhold(x))[i]);
		}

		for (j = 0; j < n_phi[i]; j++) {
			for (val = iq = 0; iq < quad[i]->n_points; iq++)
				val += phi[i][iq][j]*wdetf_qp[i][iq];
			rhs[i][j] += val;
		}
	}

#ifdef USING_MPI
	if (Rank == 1)
#endif
		if (ActualTime != t) {
			ActualTime = t;
			cout<<"Time = "<<t<<endl;
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
	FUNCNAME("main");
	
    change_msg_out(fopen("msg.info", "wt")); /* presmerovanie warningov do externeho suboru */

#ifdef USING_MPI
	if (Rank != 0)
	    change_msg_out(fopen("msg2.info", "wt")); /* presmerovanie warningov do externeho suboru */
#endif

	msg_info=2;
	ActualTime = 0;

	static int pMaskOfZorro[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, {0, 0, 1} };

	w = gsl_integration_workspace_alloc(100000);

	MyData.LoadUserDataFromFile("INIT/init.dat");
	MyData.pBuildDynamic = &build; // user has to build the global matrix (it is used at the lowest level)
	MyData.pSetUpRHS = &buildRHS; // user has to build the macro RHS, i.e., the RHS of the highest level
	MyData.pExactSolution = &u;
	MyData.pDirBoundCond = &g;
	MyData.pInitialCondition = &u;
	MyData.pMatrixMask = (int**)&(pMaskOfZorro[0][0]);

	DimOfProb = MyData.DimOfProb;

	feHMMpar MyfeHMM(&MyData);// have to be checked at the future with dimensions of the fields above !!!

	MyfeHMM.Solve();
#ifdef USING_MPI
	if (Rank == 0)
#endif
	{
		cout<<"L2Error is "<<(double)MyfeHMM.L2Error()<<endl;
		MyfeHMM.SaveSolutionToVTK("result.vt_");
		INFO(9,9,"Skoncil som\n");
	}

	gsl_integration_workspace_free(w);

	if (quad) {
		free(quad);
		free(quad_fast);
		free(phi);
		free(n_phi);
		free(wdetf_qp);
	}

#ifdef USING_MPI
	MPI_Finalize();
#endif

	return 0;
}

