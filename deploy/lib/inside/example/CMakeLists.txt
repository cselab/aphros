cmake_minimum_required(VERSION 3.3.0)
project(main C)
find_package(aphros/inside)

set(T read)
add_executable(${T} ${T}.c)
target_link_libraries(${T} inside)

set(T point)
add_executable(${T} ${T}.c)
target_link_libraries(${T} inside)

set(T point2)
add_executable(${T} ${T}.c)
target_link_libraries(${T} inside)
