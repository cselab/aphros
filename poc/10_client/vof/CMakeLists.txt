cmake_minimum_required(VERSION 3.3.0)
project(main)

find_package(aphros/hydro)
set(T imp)
add_library(${T} SHARED main.cpp)
target_link_libraries(${T} hydro)
