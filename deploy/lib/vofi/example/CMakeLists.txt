cmake_minimum_required(VERSION 3.3.0)
project(main C)
include("$ENV{CHPREFIX}/include/vofi.cmake")

set(T main)
add_executable(${T} main.c)
target_link_libraries(${T} vofi)
