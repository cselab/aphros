cmake_minimum_required(VERSION 3.3.0)
project(sph C)
include("$ENV{CHPREFIX}/include/vof.cmake")
include("$ENV{CHPREFIX}/include/vofi.cmake")
include("$ENV{CHPREFIX}/include/overlap.cmake")

add_executable(sph sph.c)
target_link_libraries(sph vof overlap)

add_executable(one one.c)
target_link_libraries(one vof)

