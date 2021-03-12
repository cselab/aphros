set(T "dumper")
aphros_add_object(T ${T} dumper.cpp)
target_link_libraries(${T} vars)

set(T "dump_hdf")
aphros_add_object(T ${T} hdf.cpp)
target_link_libraries(${T} PRIVATE use_mpi use_hdf use_dims)
if (USE_HDF)
  target_link_libraries(${T} PUBLIC hdf)
endif()

set(T "dump_raw")
aphros_add_object(T ${T} raw.cpp)
target_link_libraries(${T} PRIVATE use_mpi use_dims)

set(T "dump_xmf")
aphros_add_object(T ${T} xmf.cpp)
target_link_libraries(${T} PRIVATE parse_template use_dims use_mpi)

set(T "dump")
aphros_add_object(T ${T} dump.cpp)
target_link_libraries(${T} mesh)
