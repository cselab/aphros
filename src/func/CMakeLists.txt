set(T "primlist")
aphros_add_object(T ${T} primlist.cpp)
target_link_libraries(${T} PUBLIC format PRIVATE use_dims)

set(T "init")
aphros_add_object(T ${T} init.cpp)
target_link_libraries(${T}
    PUBLIC primlist vars suspender utilsystem
    PRIVATE use_mpi use_overlap use_dims)
if (USE_OVERLAP)
  target_link_libraries(${T} T PUBLIC overlap)
endif()

set(T "init_vel")
aphros_add_object(T ${T} init_vel.cpp)
target_link_libraries(${T} primlist vars mesh)

set(T "init_contang")
aphros_add_object(T ${T} init_contang.cpp)
target_link_libraries(${T} vars mesh)
