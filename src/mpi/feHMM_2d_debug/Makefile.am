FEHMM_DIM_OF_WORLD = 2
FEHMM_DEBUG = 1
FEHMM_MPI = 1

lib_LTLIBRARIES = libfeHMM_2d_mpi_debug.la

include $(srcdir)/../../common/Makefile_mpi.am.include

AM_CXXFLAGS = $(FEHMM_DEBUG_CXXFLAGS) $(FEHMM_CXXFLAGS)

nodist_libfeHMM_2d_mpi_debug_la_SOURCES = $(sources)
