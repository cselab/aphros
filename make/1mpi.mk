0_MPI = $(WRK)/util/subcomm.o

1_MPI = \
-D_USE_MPI_=1 \
-D_USE_HDF_=1 \

CXXFLAGS_HDF = `pkg-config --cflags hdf5-openmpi`
LDFLAGS_HDF = `pkg-config --libs hdf5-openmpi`
CXX = mpic++
CC = mpicc
