cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

set(T hook)
add_library(${T} SHARED hook.cpp)
target_link_libraries(${T} aphros)
