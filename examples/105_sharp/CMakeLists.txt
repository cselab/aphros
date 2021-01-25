cmake_minimum_required(VERSION 3.3.0)
project(main)
find_package(aphros)

set(T main)
add_executable(${T} main.cpp)
target_link_libraries(${T} aphros)

install(TARGETS ${T} DESTINATION bin)
set_target_properties(${T} PROPERTIES PREFIX "ap." OUTPUT_NAME "sharpen")

