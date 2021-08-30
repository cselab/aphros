cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

set(T hook)
add_library(${T} SHARED hook.cpp)
target_link_libraries(${T} aphros)

set(T hook)
add_library(${T} INTERFACE)
target_include_directories(${T} INTERFACE $ENV{HOME}/.local/include)
target_link_libraries(${T} INTERFACE -L$ENV{HOME}/.local/lib -llammps_serial -lmpi_stubs)

