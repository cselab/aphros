cmake_minimum_required(VERSION 3.0.0)

get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${p})

set(EXE ${p})

set(CD "../../CubismDistr")

include_directories(${CD})

add_executable(${EXE} main.cpp ${CD}/Vars.cpp)
     
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)

enable_testing()

add_test(${p} ${EXE} a b)

