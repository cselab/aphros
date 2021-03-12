set(T mesh)
aphros_add_object(T ${T} mesh.cpp)
target_link_libraries(${T} suspender use_dims use_mpi)
