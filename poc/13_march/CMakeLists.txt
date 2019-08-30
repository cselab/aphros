cmake_minimum_required(VERSION 2.8)
project(main C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(T march)
add_library(${T} SHARED
  main.c
  march.h
  )

set(CHPREFIX $ENV{CHPREFIX})

set_target_properties(${T} PROPERTIES PUBLIC_HEADER "march.h")

install(TARGETS ${T} 
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
target_include_directories(${T} PUBLIC "${CHPREFIX}/include")
target_link_libraries(${T} PUBLIC m)
install(EXPORT export FILE march.cmake DESTINATION ${CHPREFIX}/include)
