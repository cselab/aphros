cmake_minimum_required(VERSION 3.3.0)
project(main C)
include("$ENV{CHPREFIX}/include/h5.cmake")

set(T main)
add_executable(${T} main.c)
target_link_libraries(${T} h5)

set(T ba)
add_executable(${T} ba.c)
target_link_libraries(${T} h5)

set(T read)
add_executable(${T} read.c)
target_link_libraries(${T} h5)
