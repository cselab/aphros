set(T "dumper")
add_object(${T} dumper.cpp)
object_link_libraries(${T} vars)

set(T "dump_hdf")
add_object(${T} hdf.cpp)
object_link_libraries(${T} use_mpi use_hdf use_dims)
if (USE_HDF)
  object_link_libraries(${T} hdf)
endif()

set(T "dump_raw")
add_object(${T} raw.cpp)
object_link_libraries(${T} use_mpi use_dims)

set(T "dump_xmf")
add_object(${T} xmf.cpp)
object_link_libraries(${T} parse_template use_dims use_mpi)

set(T "dump_vtk")
add_object(${T} vtk.cpp)

set(T "dump")
add_object(${T} dump.cpp)
object_link_libraries(${T} mesh)
