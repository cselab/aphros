cmake_minimum_required(VERSION 3.3.0)
project(one C)
include("$ENV{CHPREFIX}/include/vof.cmake")
include("$ENV{CHPREFIX}/include/vofi.cmake")

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

set(C "main.c")
add_custom_command(
  OUTPUT "_${C}"
  COMMAND ba.env qcc -source ${C}
)

set(T main)
add_executable(${T} "_${C}")
target_link_libraries(${T} vof)
