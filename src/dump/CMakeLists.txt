set(T "dumper")
add_library(${T} dumper.cpp)
target_link_libraries(${T} vars)

set(T "dump_hdf")
add_library(${T} hdf.cpp)
target_link_libraries(${T} hdf5)
