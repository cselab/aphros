add_library(solver solver.cpp)

set(T "approx")
add_library(${T} approx.cpp)
target_link_libraries(${T} solver)

set(T "vofm")
add_library(${T} vofm.cpp)
target_link_libraries(${T} solver normal utilvof)

set(T "vof")
add_library(${T} vof.cpp)
target_link_libraries(${T} solver normal utilvof)

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
target_link_libraries(${T} solver approx)

set(T "convdiffe")
add_library(${T} convdiffe.cpp)
target_link_libraries(${T} solver approx)

set(T "convdiffe_eb")
add_library(${T} convdiffe_eb.cpp)
target_link_libraries(${T} solver approx)

set(T "convdiffv_eb")
add_library(${T} convdiffv_eb.cpp)
target_link_libraries(${T} convdiffe_eb)

set(T "convdiffvg")
add_library(${T} convdiffvg.cpp)
target_link_libraries(${T} convdiffi convdiffe)

set(T "simple")
add_library(${T} simple.cpp)
target_link_libraries(${T} convdiffvg utilconvdiff)

set(T "proj")
add_library(${T} proj.cpp)
target_link_libraries(${T} convdiffvg utilconvdiff)

set(T "proj_eb")
add_library(${T} proj_eb.cpp)
target_link_libraries(${T} convdiffv_eb utilconvdiff)
