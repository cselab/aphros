cmake_minimum_required(VERSION 3.3.0)
project(sph C)
find_package(aphros/vof)
find_package(aphros/vofi)
find_package(aphros/overlap)

add_executable(sph sph.c)
target_link_libraries(sph vof overlap)

add_executable(one one.c)
target_link_libraries(one vof)

