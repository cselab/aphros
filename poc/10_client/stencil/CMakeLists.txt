cmake_minimum_required(VERSION 3.3.0)
project(main)
include("$ENV{CHPREFIX}/include/hydro.cmake")

set(EXE main)
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} suspender)
