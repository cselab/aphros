cmake_minimum_required(VERSION 2.8)
project(vofi C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(T vofi)
add_library(${T} SHARED
                 checkconsistency.c
                 getcc.c
                 getdirs.c
                 getfh.c
                 getintersections.c
                 getlimits.c
                 getmin.c
                 getzero.c
                 integrate.c
                 interface.c
                 vofi_GL.h
                 vofi.h
                 vofi_stddecl.h)

set(CHPREFIX $ENV{CHPREFIX})

set_target_properties(${T} PROPERTIES PUBLIC_HEADER "vofi.h")

install(TARGETS ${T} 
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
target_include_directories(${T} PUBLIC "${CHPREFIX}/include")
target_link_libraries(${T} PUBLIC m)
install(EXPORT export FILE vofi.cmake DESTINATION ${CHPREFIX}/include)
