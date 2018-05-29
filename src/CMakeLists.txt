cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("mfer")

set(EXE "mfer")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  

add_subdirectory(../distr)
target_link_libraries(${EXE} distr)

add_executable(${EXE} distr/main.cpp)

set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# GCC warnings and optimizations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

# MPI
find_package(MPI REQUIRED)
add_definitions("-DMPI_ENABLE")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# C++11
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD_REQUIRED ON)

# install prefix
set(CHPREFIX $ENV{CHPREFIX})
install(TARGETS ${EXE} RUNTIME DESTINATION "${CHPREFIX}/bin")
