cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

if (APHROS_DIM4)
  set(T main)
  add_executable(${T} main.cpp)
  target_link_libraries(${T} aphros)
endif()
