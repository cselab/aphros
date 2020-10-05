cmake_minimum_required(VERSION 3.3.0)

set(CHPREFIX $ENV{CHPREFIX})

# C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

# hdf5
set(T "hdf")
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
add_library(${T} INTERFACE IMPORTED)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})

# hypre
set(T "hypreext")
set(HYPRE_DIR ${CHPREFIX})
add_library(${T} INTERFACE IMPORTED)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HYPRE_DIR}/include)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES -L${HYPRE_DIR}/lib -lHYPRE -lm)

# OpenMP
set(T "openmp")
add_library(${T} INTERFACE IMPORTED)
if (USE_OPENMP)
  find_package(OpenMP REQUIRED)
  find_package(Threads REQUIRED)
  set_property(TARGET ${T} PROPERTY 
      INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
  set_property(TARGET ${T} PROPERTY 
      INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

# Optional packages
find_package(FPZIP QUIET) # floating point compressor
