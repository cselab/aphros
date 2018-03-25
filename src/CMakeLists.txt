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
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# Source folders (add all *.cpp files)
aux_source_directory(Cubism/ SRC_CUBISM)
aux_source_directory(CubismDistr/ SRC_CUBISM_DISTR)
aux_source_directory(hydro/ SRC_HYDRO)
set(SRC ${SRC_CUBISM} ${SRC_CUBISM_DISTR} ${SRC_HYDRO})

include_directories(.)

# Add target
add_executable(${EXE} ${SRC})
     
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

# hypre
include_directories(${HYPRE_DIR}/include)
target_link_libraries(${EXE} -lHYPRE -lm)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HYPRE_DIR}/lib")

# cubism specific
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALIGNBYTES_=16") 

target_link_libraries(${EXE} ${LIBRARIES})

install(TARGETS ${EXE} RUNTIME DESTINATION $ENV{HOME}/bin)
