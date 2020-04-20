set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} overlap)

set(T "init_vel")
add_library(${T} init_vel.cpp)
