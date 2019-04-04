add_library(suspender suspender.cpp)
add_library(sysinfo sysinfo.cpp)

set(T "utilfluid")
add_library(${T} fluid.cpp)
target_link_libraries(${T} vars)
