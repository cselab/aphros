add_library(solver solver.cpp)

set(T "vof")
add_library(${T} vof.cpp)
target_link_libraries(${T} solver)
