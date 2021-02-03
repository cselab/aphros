set(T "primlist")
add_library(${T} primlist.cpp)
target_link_libraries(${T} format use_dims)

set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} primlist vars suspender use_mpi)
target_compile_definitions(${T} PUBLIC _USE_OVERLAP_=$<BOOL:${USE_OVERLAP}>)
if (USE_OVERLAP)
  target_link_libraries(${T} overlap)
endif()

set(T "init_vel")
add_library(${T} init_vel.cpp)
target_link_libraries(${T} primlist vars mesh)

set(T "init_contang")
add_library(${T} init_contang.cpp)
target_link_libraries(${T} vars mesh)

