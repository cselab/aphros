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
target_link_libraries(${T} solver)

set(T "convdiffe")
add_library(${T} convdiffe.cpp)
target_link_libraries(${T} solver)

set(T "convdiffvi")
add_library(${T} convdiffvi.cpp)
target_link_libraries(${T} convdiffi convdiffe)

set(T "simple")
add_library(${T} simple.cpp)
target_link_libraries(${T} vof tvd convdiffvi)
