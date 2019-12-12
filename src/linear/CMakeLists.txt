set(T "hypre")
add_library(${T} hypre.cpp)
target_link_libraries(${T} hypreext)

set(T "hypresub")
add_library(${T} hypresub.cpp)
target_link_libraries(${T} hypre subcomm histogram)

