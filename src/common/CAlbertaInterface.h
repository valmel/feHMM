#ifndef CALBERTA_INTERFACE_H
#define CALBERTA_INTERFACE_H

#include "alberta.h"
#include "alberta_intern.h"

#include "CUserData.h"

class CAlbertaInterface {
	MACRO_DATA *pData;
	int MeshToMacroData(MESH *);
	MESH* GetInteriorMesh();
	void InitGeometry(const char *, int); ///< nastavenie geometrie danej uzivatelom v MacroFile
	void CleanupWriteMacro(MACRO_DATA *, DOF_INT_VEC *, TRAVERSE_STACK *);
	REAL GetInnerRadius1(const EL_INFO*); ///< DIM == 1
	REAL GetInnerRadius2(const EL_INFO*); ///< DIM == 2
	REAL GetInnerRadius3(const EL_INFO*); ///< DIM == 3
	void GetNumOfBoundEl_1D();
	void GetNumOfBoundEl_2D();
	void GetNumOfBoundEl_3D();
	void SumSystemMatrices(DOF_MATRIX ***, DOF_MATRIX ***, const int **);
	void GetNumOfBoundEl();
	void FillPrimaryBoundElArray();
	void FillOppositeBoundElArray();
	void CoincidentBoundToEl(const EL_INFO*, REAL*);
	void CoincidentOppositeBoundToEl(const EL_INFO *, REAL *);
	int AreCorrespondingEls(MACRO_EL*, REAL *, int);
	void ComputeCenterOfGravity(MACRO_EL *, REAL *, int, int);
	void AssembleConstrMatrices();
	void FillOppositeDOFs(MACRO_EL *, MACRO_EL *, int);
	const DOF *pFrontDof;
public:
	// variables
	double ScaleFactor;
	int ActualLevel;
	int NumOfLevels;
	int DimOfWorld; ///< dimension of the world, i.e. if f(x,y) then DimOfWorld=2
	int DimOfProb; ///< dimension of the problem, i.e. if f=(f1,f2,f3) then DimOfProb=3
	MESH* pMesh; ///< mesh; initialized in InitReference or InitGeometry
	MESH* pRestrictedMesh; ///< mesh restricted to interior elements
	int NumOfRefin; ///< number of refinements
	int QuadDeg; ///< degree of quadrature
	const QUAD* pQuad; ///< quadrature
	int *LagDeg; ///< degree of finite elements, Lagrange only!!!
	const BAS_FCTS **Lagrange;
	const FE_SPACE **FEspace; ///< space of finite elements
	const FE_SPACE *ProblemFEspace; ///< space of finite elements
	const FE_SPACE **RestrictedFEspace; ///< space of finite elements on restricted mesh
	const FE_SPACE *LinearFEspace; ///< space of finite elements used on the boundary (periodic case)
	int *NumOfBasFct; ///< number of basis functions
	int *IdxOfDofs; ///< number of dofs in every FEspace
	DOF **Dofs;
	DOF **RestrictedDofs;
	DOF_MATRIX ***Matrix;
	DOF_MATRIX ***RestrictedMatrix;
	DOF_MATRIX ***StaticMatrix;
	DOF_MATRIX ***StaticRestrictedMatrix;
	DOF_MATRIX *ProblemMatrix;
	DOF_REAL_VEC *ProblemSol;
	DOF_REAL_VEC *ProblemRHS;
	DOF_REAL_VEC ***ProblemUh;
	DOF_REAL_VEC **Uh;
	DOF_REAL_VEC **Fh;
	DOF_INT_VEC *pOppositeDOFs;
	REAL **LocalFh;
	DOF_REAL_VEC **pRestrictedVertexIndex; ///< mapping from the restricted matrix to extended
	DOF_REAL_VEC **pExtendedVertexIndex; ///< mapping from extended matrix to restricted (-1 if no mapping exists)
	REAL ****LocalMatrix;
	DOF_MATRIX **ConstrMatrix;
	DOF_MATRIX **ConstrMatrixT;

	REAL_D x, y;
	const DOF *(**GetGlobDof)(const EL *, const DOF_ADMIN *, DOF *);
	const S_CHAR *(**GetBound)(const EL_INFO *, S_CHAR *);
	REAL *pVertices; ///< world coordinates of reference elements ->To AuxiliaryData...
	REAL *pRestrictedVertices; ///< world coordinates of reference elements ->To AuxiliaryData...
	int NumOfBoundEl;
	MACRO_EL**** pBoundaryEl;
	MACRO_EL** pAllBoundaryEl;
	int BoundCond;

	// functions
	virtual void Initialize(CUserData *, int, int);
	virtual void Free();
	void TransformMeshes(REAL* dx); ///< transformacia mriezky
	int IsInside(const EL_INFO* el_info);
	int GetProblemDof(int, int);
	void ClearLocalMatrix();
	void ClearMatrix();
	void ClearRestrictedMatrix();
	int FillRestrictedVertexIndex();
	void SumMatrices(const int**);
	void SumRestrictedMatrices(const int **);
	REAL ComputeMeshVolume(MESH*);
	void AssembleProblemMatrix(const int **, const DOF_MATRIX ***); ///<gues what
	void ExtendRestrictedMatrix(int, int, const DOF_MATRIX ***);
	void PostProcessMatrix(REAL);
	void CreateStaticMatrices();
	REAL GetInnerRadius(const EL_INFO*); ///< el_info with coordinates
	CAlbertaInterface() :
		pData(NULL) {
	}
	;
	virtual ~CAlbertaInterface() {
	}
	;
	MESH* GetRestrictedMesh(REAL*, REAL*, MESH*);
	MACRO_DATA* GetReferenceMacroData(int Dim); ///< function returns reference qube depending on Dim
};

#endif // #ifdef CALBERTA_INTERFACE_H
