cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("main")

set(EXE "main")

include_directories(.)

# Add target
add_executable(${EXE} main.cpp reconst.cpp)
