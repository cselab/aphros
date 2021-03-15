0_MPI = $(WRK)/util/subcomm_dummy.o

1_MPI = \
-D_USE_MPI_=0 \
-D_USE_HDF_=0 \

CXXFLAGS_HDF =
LDFLAGS_HDF =
CXX = g++
CC = c99
