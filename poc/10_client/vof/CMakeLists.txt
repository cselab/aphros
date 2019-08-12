cmake_minimum_required(VERSION 3.3.0)
project("main")

include("$ENV{CHPREFIX}/include/hydro.cmake")
set(T "imp")
add_library(${T} main.cpp)
target_link_libraries(${T} hydro)

set(EXE main)
add_executable(${EXE} main.c)
target_link_libraries(${EXE} imp)
