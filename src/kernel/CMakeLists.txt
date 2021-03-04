set(T "hydro")
add_library(${T} hydro.cpp)
target_link_libraries(${T} aphros use_dims)
