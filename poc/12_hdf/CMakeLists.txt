cmake_minimum_required(VERSION 3.3.0)
project("hdf" C)

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

set(T hdf)
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)

if(NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "no parallel HDF5")
endif()

add_library(${T} INTERFACE IMPORTED)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})

set(T main)
add_executable(${T} main.c)
target_link_libraries(${T} hdf)

set(T ehdf)
add_executable(${T} hdf.c)
target_link_libraries(${T} hdf)
