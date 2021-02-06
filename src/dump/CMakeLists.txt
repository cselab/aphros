set(T "dumper")
add_library(${T} dumper.cpp)
target_link_libraries(${T} vars)

set(T "dump_hdf")
add_library(${T} hdf.cpp)
target_link_libraries(${T} PRIVATE use_mpi use_hdf use_dims)
if (USE_HDF)
  target_link_libraries(${T} PUBLIC hdf)
endif()

set(T "dump_raw")
add_library(${T} raw.cpp)
target_link_libraries(${T} PRIVATE use_mpi use_dims)

set(T "dump")
add_library(${T} dump.cpp)
target_link_libraries(${T} mesh)
