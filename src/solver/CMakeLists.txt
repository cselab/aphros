set(T "solver")
add_library(${T} solver.cpp)
target_link_libraries(${T} timer)

set(T "approx")
add_library(${T} approx.cpp)
target_link_libraries(${T} solver)

set(T "embed")
add_library(${T} embed.cpp)
target_link_libraries(${T} suspender)

set(T "approx_eb")
add_library(${T} approx_eb.cpp)
target_link_libraries(${T} approx embed vars inside)

set(T "vofm")
add_library(${T} vofm.cpp)
target_link_libraries(${T} solver normal utilvof)

set(T "vof")
add_library(${T} vof.cpp)
target_link_libraries(${T} solver normal utilvof approx_eb)

set(T "tvd")
add_library(${T} tvd.cpp)
target_link_libraries(${T} solver)

set(T "normal")
add_library(${T} normal.cpp)
target_link_libraries(${T} approx)

set(T "curv")
add_library(${T} curv.cpp)
target_link_libraries(${T} partstrmeshm normal vof vofm)

set(T "partstrmeshm")
add_library(${T} partstrmeshm.cpp)
target_link_libraries(${T} dumper)

set(T "convdiffi")
add_library(${T} convdiffi.cpp)
target_link_libraries(${T} solver approx approx_eb)

set(T "convdiffe")
add_library(${T} convdiffe.cpp)
target_link_libraries(${T} solver approx approx_eb)

set(T "convdiffvg")
add_library(${T} convdiffvg.cpp)
target_link_libraries(${T} convdiffi convdiffe)

set(T "simple")
add_library(${T} simple.cpp)
target_link_libraries(${T} convdiffvg utilconvdiff)

set(T "proj")
add_library(${T} proj.cpp)
target_link_libraries(${T} utilconvdiff)

set(T "proj_eb")
add_library(${T} proj_eb.cpp)
target_link_libraries(${T} utilconvdiff)
