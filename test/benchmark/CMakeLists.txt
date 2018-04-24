cmake_minimum_required(VERSION 3.0.0)

get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${p})

set(CMAKE_BUILD_TYPE release)

set(EXE t.${p})

set(S "../../src")

include_directories(${S})

add_executable(${EXE} main.cpp ${S}/util/suspender.cpp timer.cpp)
     
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)

enable_testing()

add_test(${p} ${EXE})

