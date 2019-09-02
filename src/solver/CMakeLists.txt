add_library(solver solver.cpp)

set(T "vofm")
add_library(${T} vofm.cpp)
target_link_libraries(${T} solver normal partstrmeshm march utilvof)

set(T "vof")
add_library(${T} vof.cpp)
target_link_libraries(${T} solver normal partstrmesh utilvof)

set(T "tvd")
add_library(${T} tvd.cpp)
target_link_libraries(${T} solver)

set(T "normal")
add_library(${T} normal.cpp)

set(T "partstrmesh")
add_library(${T} partstrmesh.cpp)
target_link_libraries(${T} dumper)

set(T "partstrmeshm")
add_library(${T} partstrmeshm.cpp)
target_link_libraries(${T} dumper)

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
target_link_libraries(${T} vof vofm tvd convdiffvi)

set(T "multi")
add_library(${T} multi.cpp)
