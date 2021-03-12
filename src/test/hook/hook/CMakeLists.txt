cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

set(T hook)
aphros_add_object(T ${T} SHARED hook.cpp)
target_link_libraries(${T} aphros)

