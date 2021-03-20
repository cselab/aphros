set(T hydro_post)
add_object(${T} hydro_post.cpp)
object_link_libraries(${T} PUBLIC dump visual PRIVATE use_dims use_mpi)

set(T hydro)
add_object(${T} hydro.cpp)
object_link_libraries(${T}
  PUBLIC convdiffe convdiffi distrsolver utildistr proj simple dump hydro_post
  PRIVATE use_dims use_mpi)

