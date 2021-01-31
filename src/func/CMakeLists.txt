set(T "primlist")
add_library(${T} primlist.cpp)
target_link_libraries(${T} format use_dims)

set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} primlist overlap vars suspender)

set(T "init_vel")
add_library(${T} init_vel.cpp)
target_link_libraries(${T} primlist vars)

set(T "init_contang")
add_library(${T} init_contang.cpp)
target_link_libraries(${T} vars)

