cmake_minimum_required(VERSION 3.3.0)
project(vof C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CHPREFIX $ENV{CHPREFIX})
set(CMAKEPREFIX ${CHPREFIX}/cmake/aphros)

find_package("aphros/vofi")
set(T vof)
add_library(${T} SHARED vof.h vof.c)
target_link_libraries(${T} PUBLIC vofi)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "vof.h")
install(TARGETS ${T}
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
install(EXPORT export FILE vofConfig.cmake DESTINATION ${CMAKEPREFIX})
