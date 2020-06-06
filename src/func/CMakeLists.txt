set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} overlap vars suspender)

set(T "init_vel")
add_library(${T} init_vel.cpp)
target_link_libraries(${T} vars)

set(T "init_contang")
add_library(${T} init_contang.cpp)
target_link_libraries(${T} vars)
