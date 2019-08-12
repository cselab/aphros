set(T "hypre")

add_library(${T} hypre.cpp)
target_link_libraries(${T} hypreext)

