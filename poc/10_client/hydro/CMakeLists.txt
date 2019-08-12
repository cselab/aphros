cmake_minimum_required(VERSION 3.3.0)

project("main")

set(EXE "main")

include("$ENV{CHPREFIX}/include/client.cmake")
include("$ENV{CHPREFIX}/include/hydro.cmake")

# binary
add_executable(${EXE} main.cpp)
include_directories(${EXE} ${CHPREFIX}/include/hydro)
link_directories(${EXE} ${CHPREFIX}/lib)
target_link_libraries(${EXE} hydro)
