#ifndef ERROR_RETURN_CODES_H
#define ERROR_RETURN_CODES_H

/// the codes of the return values of functions
enum EErrorReturnCodes
{
	ERC_OK=0, ///< no problem
	ERC_NUM_OF_LEVELS_MISSING, ///<
	ERC_DIM_OF_WORLD_MISSING, ///<
	ERC_DIM_OF_PROB_MISSING, ///<
	ERC_REFINEMENT_MISSING, ///<
	ERC_QUADRATURE_MISSING, ///<
	ERC_SCALE_FACTOR_MISSING, ///<
	ERC_FE_DEG_MISSING, ///<
	ERC_BUILD_FUNC_MISSING, ///<
	ERC_SETUPRHS_FUNC_MISSING, ///<
	ERC_MATRIX_NOT_ASSEMBLED, ///< solve used without a prior assemblage
	ERC_ERROR_OPEN_VTK_FILE, ///< urcite mame definovany vystupny vtk subor?
	ERC_ERROR_MACRO_FILE_MISSING, ///< Alberta's macro geometry file not given
	ERC_SOLVER_NOT_SET, ///< user probably forgot to set solver
	ERC_IS_BI_LIN_FORM_DEPENDENT_ON_PREV_STEP_MISSING, ///< we are solving a time problem and do not know if the solution from the previous time step is needed to assemble
	ERC_TIME_INFO_NOT_COMPLETE, ///< user probably forgot to set StartTime, EndTime, TimeStep or NumTimeStep
	ERC_ZERO_NUM_OF_TIME_STEP, ///< number of steps is zero
	ERC_MATRIX_MASK_NOT_SET, ///< matrix mask not set
	ERC_LEVELS_NOT_INITIALIZED, ///< this should not happend if yes, problem in the library
	ERC_SOLUTION_MISSING, ///< in save vtk is solution missing
};

#endif  // ifndef ERROR_RETURN_CODES_H
