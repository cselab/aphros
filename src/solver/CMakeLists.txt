add_library(solver solver.cpp)

set(T "vof")
add_library(${T} vof.cpp)
target_link_libraries(${T} solver normal)

set(T "tvd")
add_library(${T} tvd.cpp)
target_link_libraries(${T} solver)

set(T "normal")
add_library(${T} normal.cpp)

set(T "convdiffi")
add_library(${T} convdiffi.cpp)

set(T "convdiffvi")
add_library(${T} convdiffvi.cpp)
target_link_libraries(${T} convdiffi)
