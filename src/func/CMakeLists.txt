set(T "primlist")
add_library(${T} primlist.cpp)
target_link_libraries(${T} PUBLIC format PRIVATE use_dims)

set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} PUBLIC primlist vars suspender utilsystem PRIVATE use_mpi
    use_overlap use_dims)
if (USE_OVERLAP)
  target_link_libraries(${T} PUBLIC overlap)
endif()

set(T "init_vel")
add_library(${T} init_vel.cpp)
target_link_libraries(${T} primlist vars mesh)

set(T "init_contang")
add_library(${T} init_contang.cpp)
target_link_libraries(${T} vars mesh)

