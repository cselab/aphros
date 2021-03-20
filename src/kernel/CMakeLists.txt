set(T hydro)
add_object(${T} hydro.cpp)
object_link_libraries(${T}
  PUBLIC convdiffe convdiffi distrsolver utildistr proj simple dump hydro_post
  PRIVATE use_dims use_mpi)

