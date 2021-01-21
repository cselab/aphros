cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

set(EXE main)
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} aphros)

