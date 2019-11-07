cmake_minimum_required(VERSION 3.3.0)
project("hdf" C)

#find_package(MPI REQUIRED)

#set(HDF5_PREFER_PARALLEL on)
#find_package(HDF5 REQUIRED COMPONENTS C HL)
#set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

#if(NOT HDF5_IS_PARALLEL)
#  message(FATAL_ERROR "no parallel HDF5")
#endif()

#message(1 "${HDF5_FOUND}")
#message(2 "${HDF5_VERSION}")
#message(3 "${HDF5_INCLUDE_DIRS}")
#message(4 "${HDF5_INCLUDE_DIR}")
#message(5 "${HDF5_DEFINITIONS}")
#message(6 "${HDF5_LIBRARIES}")
#message(7 "${HDF5_HL_LIBRARIES}")

set(T main)
#target_link_libraries(${T} PUBLIC hdf)
add_executable(${T} main.c)
