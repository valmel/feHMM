FEHMM_DIM_OF_WORLD = 3
FEHMM_DEBUG = 0
FEHMM_MPI = 1

lib_LTLIBRARIES = libfeHMM_3d_mpi.la

include $(srcdir)/../../common/Makefile_mpi.am.include

AM_CXXFLAGS = $(FEHMM_RELEASE_CXXFLAGS) $(FEHMM_CXXFLAGS)

nodist_libfeHMM_3d_mpi_la_SOURCES = $(sources)
