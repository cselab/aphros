set(T "primlist")
add_object(${T} primlist.cpp)
object_link_libraries(${T} format use_dims)

set(T "init")
add_object(${T} init.cpp)
object_link_libraries(${T}
  primlist
  suspender
  use_mpi use_dims
  utilsystem
  vars
)
object_link_libraries(${T} overlap)

set(T "init_vel")
add_object(${T} init_vel.cpp)
object_link_libraries(${T} primlist vars mesh)

set(T "init_contang")
add_object(${T} init_contang.cpp)
object_link_libraries(${T} vars mesh)
