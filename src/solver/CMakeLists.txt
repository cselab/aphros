set(T "solver")
add_object(${T} solver.cpp)
object_link_libraries(${T} timer)

set(T "approx")
add_object(${T} approx.cpp)
object_link_libraries(${T} solver suspender use_dims use_mpi)

set(T "embed")
add_object(${T} embed.cpp)
object_link_libraries(${T} suspender mesh dump_vtk)

set(T "approx_eb")
add_object(${T} approx_eb.cpp)
object_link_libraries(${T} approx embed vars  use_dims primlist inside)

set(T "vofm")
add_object(${T} vofm.cpp)
object_link_libraries(${T} solver normal utilvof approx_eb dump_hdf)

set(T "vof")
add_object(${T} vof.cpp)
object_link_libraries(${T} solver normal utilvof approx_eb use_dims)

set(T "tracer")
add_object(${T} tracer.cpp)
object_link_libraries(${T} approx_eb use_dims use_mpi)

set(T "particles")
add_object(${T} particles.cpp)
object_link_libraries(${T} mesh use_dims)

set(T "normal")
add_object(${T} normal.cpp)
object_link_libraries(${T} approx_eb use_dims use_avx use_mpi)

set(T "curv")
add_object(${T} curv.cpp)
object_link_libraries(${T} partstrmeshm normal vof vofm)

set(T "partstrmeshm")
add_object(${T} partstrmeshm.cpp)
object_link_libraries(${T} dumper suspender use_mpi use_dims dump_vtk)

set(T "convdiffi")
add_object(${T} convdiffi.cpp)
object_link_libraries(${T} solver approx approx_eb)

set(T "convdiffe")
add_object(${T} convdiffe.cpp)
object_link_libraries(${T} solver approx approx_eb)

set(T "convdiffvg")
add_object(${T} convdiffvg.cpp)
object_link_libraries(${T} convdiffi convdiffe)

set(T "simple")
add_object(${T} simple.cpp)
object_link_libraries(${T} convdiffvg utilconvdiff)

set(T "proj")
add_object(${T} proj.cpp)
object_link_libraries(${T} utilconvdiff use_dims)

set(T "proj_eb")
add_object(${T} proj_eb.cpp)
object_link_libraries(${T} utilconvdiff use_dims)

set(T "fluid_dummy")
add_object(${T} fluid_dummy.cpp)
object_link_libraries(${T} vars suspender posthook_default approx_eb)

set(T "electro")
add_object(${T} electro.cpp)
object_link_libraries(${T} approx_eb use_dims use_mpi)

