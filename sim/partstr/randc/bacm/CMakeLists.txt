cmake_minimum_required(VERSION 3.3.0)
set(CHPREFIX $ENV{CHPREFIX})
project(one C)
include("$ENV{CHPREFIX}/include/vof.cmake")
include("$ENV{CHPREFIX}/include/vofi.cmake")

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

set(C "main.c")
add_custom_target(qcc COMMAND ba.env qcc -source ${C} BYPRODUCTS _${C})
add_custom_target(u COMMAND ba.env ba.u
  "${CHPREFIX}/include/curv/curv.conf"
  "${CHPREFIX}/include/bashape.conf")

set(T main.s)
add_executable(${T} "_${C}")
add_dependencies(qcc u)
add_dependencies(${T} qcc)
target_link_libraries(${T} vof)
