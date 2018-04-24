cmake_minimum_required(VERSION 3.0.0)

get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${p})

set(EXE t.${p})

set(H "../../src/hydro")
set(CD "../../src/CubismDistr")
set(S "../../src")

include_directories(${S})

add_executable(${EXE} main.cpp 
    ${H}/suspender.cpp 
    ${H}/hypre.cpp 
    ${H}/dumper.cpp
    ${CD}/Vars.cpp
    ${CD}/Local.cpp
    ${CD}/Cubism.cpp
    ${CD}/Interp.cpp
    ${CD}/DistrSolver.cpp
    )
     
# OpenMP 
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS 
    "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# HDF5
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
add_definitions(" -D_USE_HDF_")
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(${EXE} ${HDF5_LIBRARIES})

# hypre
include_directories(${HYPRE_DIR}/include)
target_link_libraries(${EXE} -lHYPRE -lm)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HYPRE_DIR}/lib")

# cubism specific
set(BS 16)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_BLOCKSIZE_=${BS} -D_BLOCKSIZEX_=${BS} -D_BLOCKSIZEY_=${BS} -D_BLOCKSIZEZ_=${BS}") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALIGNBYTES_=16") 


set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)

enable_testing()

add_test(${p} ${EXE})
