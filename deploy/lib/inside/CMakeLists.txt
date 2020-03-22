cmake_minimum_required(VERSION 2.8)
project(main C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(T inside)
add_library(${T} SHARED
  main.c
  march.h
  )

set(CHPREFIX $ENV{CHPREFIX})
set(CMAKEPREFIX ${CHPREFIX}/cmake/aphros)

set_target_properties(${T} PROPERTIES PUBLIC_HEADER "inside.h")

install(TARGETS ${T} 
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
target_include_directories(${T} PUBLIC "${CHPREFIX}/include")
target_link_libraries(${T} PUBLIC m)
install(EXPORT export FILE marchConfig.cmake DESTINATION ${CMAKEPREFIX})
