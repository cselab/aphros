cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros/hydro)

set(T posthook)
add_library(${T} SHARED posthook.cpp)
target_link_libraries(${T} suspender vars)

