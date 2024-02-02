# AC_PREREQ(2.61)
AC_INIT([feHMM], [1.0], [bug@fehmm.org])
AC_CONFIG_SRCDIR([src/common/feHMMpar.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LIBTOOL

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# Checks if alberta paths are set. They have to be!

if { test ! -n "${ALBERTA_INCLUDE_PATH+1}"; } || { test ! -n "${ALBERTA_DATA_PATH+1}"; } || { test ! -n "${ALBERTA_LIB_PATH+1}"; }; then
AC_MSG_NOTICE([Please, define Alberta's paths (in .bashrc for example)])
exit
fi

#OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ALBERTA_LIB_PATH}
AC_CHECK_LIB([alberta_1d],[el_det_1d])
AC_CHECK_LIB([alberta_2d],[el_det_2d])
AC_CHECK_LIB([alberta_3d],[el_det_3d])
#LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}

m4_include(m4/acx_mpi.m4)
AC_ARG_WITH([mpi],[AS_HELP_STRING([--with-mpi],[if you want mpi libraries to be build])],
            [with_mpi=yes],
            [with_mpi=no])

# Checks mpi consistency if --with-mpi is present 

if test "$with_mpi" = yes; then
AC_LANG_CPLUSPLUS
ACX_MPI([AC_MSG_NOTICE([mpi libraries will be build])],[with_mpi=no])
AC_LANG_C
ACX_MPI([],[])
fi

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow sqrt])

FFLAGS=""
CFLAGS=""
CXXFLAGS=""

FEHMM_CXXFLAGS="-I${ALBERTA_INCLUDE_PATH} -I${ALBERTA_INCLUDE_PATH}/../../alberta/src/Common -I${ALBERTA_INCLUDE_PATH}/../../alberta/src/1d -I${ALBERTA_INCLUDE_PATH}/../../alberta/src/2d -I${ALBERTA_INCLUDE_PATH}/../../alberta/src/3d"
AC_SUBST(FEHMM_CXXFLAGS)

# Here standard flags are defined. They could be overloaded.

if test -z "${FEHMM_RELEASE_CXXFLAGS}"; then
	if test "$CXX"="g++"; then
		FEHMM_RELEASE_CXXFLAGS="-O3"
	else
		FEHMM_RELEASE_CXXFLAGS="-O"
	fi
fi

AC_ARG_VAR([FEHMM_RELEASE_CXXFLAGS],[Compiler-flags used to create the optimised libraries])

if test -z "${FEHMM_DEBUG_CXXFLAGS}"; then
	if test "$CXX"="g++"; then
		FEHMM_DEBUG_CXXFLAGS="-ggdb3"
       	else
		FEHMM_DEBUG_CXXFLAGS="-g"
	fi
fi

AC_ARG_VAR([FEHMM_DEBUG_CXXFLAGS],[Compiler-flags used to create the debug-enabled libraries])

AM_CONDITIONAL([COND_MPI], [test "$with_mpi" = yes])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/common/Makefile
  src/feHMM_1d/Makefile
  src/feHMM_1d_debug/Makefile
  src/feHMM_2d/Makefile
  src/feHMM_2d_debug/Makefile
  src/feHMM_3d/Makefile
  src/feHMM_3d_debug/Makefile
  src/mpi/Makefile
  src/mpi/feHMM_1d/Makefile
  src/mpi/feHMM_1d_debug/Makefile
  src/mpi/feHMM_2d/Makefile
  src/mpi/feHMM_2d_debug/Makefile
  src/mpi/feHMM_3d/Makefile
  src/mpi/feHMM_3d_debug/Makefile
  demo/elliptic/Makefile
])
AC_OUTPUT