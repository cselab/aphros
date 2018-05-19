cmake_minimum_required(VERSION 3.0.0)

get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${p})

set(EXE t.${p})

set(S "../../src")

include_directories(${S})

add_executable(${EXE} main.cpp ${S}/linear/hypre.cpp)
     
# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# install prefix
set(CHPREFIX $ENV{CHPREFIX})

# hypre
set(HYPRE_DIR ${CHPREFIX})
include_directories(${HYPRE_DIR}/include)
target_link_libraries(${EXE} -lHYPRE -lm)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HYPRE_DIR}/lib")

set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)

enable_testing()

add_test(${p} ${EXE})
