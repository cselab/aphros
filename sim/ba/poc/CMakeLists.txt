cmake_minimum_required(VERSION 3.3.0)
project(one C)
include("$ENV{CHPREFIX}/include/vof.cmake")
include("$ENV{CHPREFIX}/include/vofi.cmake")

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

# git revision
set(P ".")
set(C "${P}/_hello.c")
add_custom_command(
  OUTPUT "${C}"
  COMMAND qcc -source hello.c
  WORKING_DIRECTORY "${P}"
)

set(T hello)
add_executable(${T} "${C}")
target_link_libraries(${T} vof)
