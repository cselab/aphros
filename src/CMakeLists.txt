cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("mfer")

set(EXE "mfer")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  

# OpenMP support                              
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              

# GCC warnings and optimizations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

set(S ".")

include_directories(${S})

add_executable(${EXE} 
    ${S}/distr/main.cpp
    ${S}/util/suspender.cpp 
    ${S}/dump/dumper.cpp
    ${S}/parse/vars.cpp
    ${S}/distr/local.cpp
    ${S}/distr/cubism.cpp
    ${S}/parse/interp.cpp
    ${S}/distr/distrsolver.cpp
    )

set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# C++11
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD_REQUIRED ON)

# MPI
find_package(MPI REQUIRED)
add_definitions("-DMPI_ENABLE")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# HDF5
if (1) 
  set(HDF5_PREFER_PARALLEL on)
  find_package(HDF5 REQUIRED)
  add_definitions(" -D_USE_HDF_")
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(${EXE} ${HDF5_LIBRARIES})
endif()

# install prefix
set(CHPREFIX $ENV{CHPREFIX})

# cubism specific
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALIGNBYTES_=16") 

add_subdirectory(linear)
target_link_libraries(${EXE} hypre)


install(TARGETS ${EXE} RUNTIME DESTINATION "${CHPREFIX}/bin")
