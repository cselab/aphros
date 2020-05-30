cmake_minimum_required(VERSION 3.3.0)
project(main C)
find_package(aphros/h5)

set(T main)
add_executable(${T} main.c)
target_link_libraries(${T} h5)

set(T ba)
add_executable(${T} ba.c)
target_link_libraries(${T} h5)

set(T dump)
add_executable(${T} dump.c)
target_link_libraries(${T} h5)

set(T write_serial)
add_executable(${T} write_serial.c)
target_link_libraries(${T} h5)
