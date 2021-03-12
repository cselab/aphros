set(T "hydro")
aphros_add_object(T ${T} hydro.cpp)
target_link_libraries(${T} aphros use_dims use_mpi)
