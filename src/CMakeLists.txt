cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("mfer")

set(EXE "mfer")

add_executable(${EXE} main.cpp)

set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# install prefix
set(CHPREFIX $ENV{CHPREFIX})
install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")

# MPI
find_package(MPI REQUIRED)
add_definitions("-DMPI_ENABLE")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# GCC warnings and optimizations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

# C++11
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD_REQUIRED ON)

# OpenMP 
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              

include_directories(.)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(parse)
add_subdirectory(dump)

target_link_libraries(${EXE} distr)

add_subdirectory(test)
enable_testing()
