set(T mesh)
add_object(${T} mesh.cpp)
object_link_libraries(${T} suspender use_dims use_mpi)
