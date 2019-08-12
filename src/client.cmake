cmake_minimum_required(VERSION 3.3.0)

set(CHPREFIX $ENV{CHPREFIX})

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wno-deprecated-copy")

# hdf5
set(T "hdf")
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
add_library(${T} INTERFACE IMPORTED)
target_include_directories(${T} PUBLIC INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${T} INTERFACE ${HDF5_LIBRARIES})

# hypre
set(T "hypreext")
set(HYPRE_DIR ${CHPREFIX})
add_library(${T} INTERFACE IMPORTED)
target_include_directories(${T} PUBLIC INTERFACE ${HYPRE_DIR}/include)
target_link_libraries(${T} INTERFACE -L${HYPRE_DIR}/lib -lHYPRE -lm)
