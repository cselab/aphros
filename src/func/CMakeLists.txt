set(T "init")
add_library(${T} init.cpp)
target_link_libraries(${T} overlap)
