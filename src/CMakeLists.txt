cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("mfer")

set(EXE "mfer")

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# GCC warnings and optimizations
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(${EXE} main.cpp)
set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# install prefix
set(CHPREFIX $ENV{CHPREFIX})
install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")

include_directories(.)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(parse)
add_subdirectory(dump)

target_link_libraries(${EXE} distr)

add_subdirectory(test)
enable_testing()
