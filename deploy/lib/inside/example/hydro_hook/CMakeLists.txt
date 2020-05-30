cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros/hydro)
find_package(aphros/inside)

set(T hook)
add_library(${T} SHARED hook.cpp)
target_link_libraries(${T} vars inside)

