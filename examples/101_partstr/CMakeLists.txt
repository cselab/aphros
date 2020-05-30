cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros/hydro)

set(EXE main)
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} distrbasic curv)

