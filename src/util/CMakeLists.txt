add_library(suspender suspender.cpp)
add_library(sysinfo sysinfo.cpp)

set(T "utilhydro")
add_library(${T} hydro.cpp)
target_link_libraries(${T} vars)

set(T "utilconvdiff")
add_library(${T} convdiff.cpp)
target_link_libraries(${T} convdiffvg)

set(T "utilvof")
add_library(${T} vof.cpp)
target_link_libraries(${T} march)

set(T "events")
add_library(${T} events.cpp)
target_link_libraries(${T} vars parser)

set(T "timer")
add_library(${T} timer.cpp)

add_library(gitgen gitgen.cpp)
add_dependencies(gitgen gitrev)

set(T "git")
add_library(${T} git.cpp)
target_link_libraries(${T} gitgen)
