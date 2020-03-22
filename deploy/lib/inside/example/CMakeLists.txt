cmake_minimum_required(VERSION 3.3.0)
project(main C)
find_package(aphros/inside)

set(T read)
add_executable(${T} read.c)
target_link_libraries(${T} march)
