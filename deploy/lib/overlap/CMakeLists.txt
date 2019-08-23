cmake_minimum_required(VERSION 2.8)
project(overlap CXX)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(T overlap)
add_library(${T} SHARED
  overlap.cpp
  overlap.h
  overlap.hpp
)
set(CHPREFIX $ENV{CHPREFIX})
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "overlap.h")
install(TARGETS ${T} 
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
target_include_directories(${T} PUBLIC "${CHPREFIX}/include")
message(cmake eigen3 $ENV{EIGEN3_INCLUDE_DIR})
target_include_directories(${T} PRIVATE $ENV{EIGEN3_INCLUDE_DIR})
set_property(TARGET ${T} PROPERTY CXX_STANDARD 11)
target_link_libraries(${T} PUBLIC m)
target_link_libraries(${T} PUBLIC stdc++)
install(EXPORT export FILE overlap.cmake DESTINATION ${CHPREFIX}/include)
