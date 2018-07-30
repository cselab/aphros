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

# warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Werror")

add_executable(${EXE} main.cpp)
set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# install prefix
set(CHPREFIX $ENV{CHPREFIX})
install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")

include_directories(.)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(solver)
add_subdirectory(parse)
add_subdirectory(dump)

target_link_libraries(${EXE} distr solver vof)

add_subdirectory(test)
enable_testing()
