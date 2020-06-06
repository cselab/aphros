set(T "dumper")

add_library(${T} dumper.cpp)
target_link_libraries(${T} vars)
