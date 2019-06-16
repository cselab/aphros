add_library(suspender suspender.cpp)
add_library(sysinfo sysinfo.cpp)

set(T "utilfluid")
add_library(${T} fluid.cpp)
target_link_libraries(${T} vars)

set(T "events")
add_library(${T} events.cpp)
target_link_libraries(${T} vars parser)

set(T "timer")
add_library(${T} timer.cpp)

add_library(gitgen gitgen.cpp)

set(T "git")
add_library(${T} git.cpp)
target_link_libraries(${T} gitgen)
