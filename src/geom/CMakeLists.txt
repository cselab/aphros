set(T mesh)
add_library(${T} mesh.cpp)
target_link_libraries(${T} suspender use_dims use_mpi)
