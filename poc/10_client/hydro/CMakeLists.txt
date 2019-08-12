cmake_minimum_required(VERSION 3.3.0)

project("main")

set(EXE "main")

set(CHPREFIX $ENV{CHPREFIX})

include("${CHPREFIX}/include/client.cmake")

# binary
add_executable(${EXE} main.cpp)
include_directories(${EXE} ${CHPREFIX}/include/hydro)
link_directories(${EXE} ${CHPREFIX}/lib)
target_link_libraries(${EXE} hydro)
