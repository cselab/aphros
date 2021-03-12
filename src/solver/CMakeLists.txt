set(T "solver")
aphros_add_object(T ${T} solver.cpp)
target_link_libraries(${T} timer)

set(T "approx")
aphros_add_object(T ${T} approx.cpp)
target_link_libraries(${T} solver suspender use_dims use_mpi)

set(T "embed")
aphros_add_object(T ${T} embed.cpp)
target_link_libraries(${T} suspender mesh)

set(T "approx_eb")
aphros_add_object(T ${T} approx_eb.cpp)
target_link_libraries(${T} PRIVATE approx embed vars  use_dims primlist PUBLIC inside)

set(T "vofm")
aphros_add_object(T ${T} vofm.cpp)
target_link_libraries(${T} solver normal utilvof approx_eb dump_hdf)

set(T "vof")
aphros_add_object(T ${T} vof.cpp)
target_link_libraries(${T} solver normal utilvof approx_eb use_dims)

set(T "tracer")
aphros_add_object(T ${T} tracer.cpp)
target_link_libraries(${T} PUBLIC approx_eb PRIVATE use_dims use_mpi)

set(T "particles")
aphros_add_object(T ${T} particles.cpp)
target_link_libraries(${T} PUBLIC mesh PRIVATE use_dims)

set(T "normal")
aphros_add_object(T ${T} normal.cpp)
target_link_libraries(${T} PUBLIC approx_eb use_dims PRIVATE use_avx use_mpi)

set(T "curv")
aphros_add_object(T ${T} curv.cpp)
target_link_libraries(${T} partstrmeshm normal vof vofm)

set(T "partstrmeshm")
aphros_add_object(T ${T} partstrmeshm.cpp)
target_link_libraries(${T} PUBLIC dumper suspender PRIVATE use_mpi use_dims)

set(T "convdiffi")
aphros_add_object(T ${T} convdiffi.cpp)
target_link_libraries(${T} solver approx approx_eb)

set(T "convdiffe")
aphros_add_object(T ${T} convdiffe.cpp)
target_link_libraries(${T} solver approx approx_eb)

set(T "convdiffvg")
aphros_add_object(T ${T} convdiffvg.cpp)
target_link_libraries(${T} convdiffi convdiffe)

set(T "simple")
aphros_add_object(T ${T} simple.cpp)
target_link_libraries(${T} convdiffvg utilconvdiff)

set(T "proj")
aphros_add_object(T ${T} proj.cpp)
target_link_libraries(${T} PRIVATE utilconvdiff PUBLIC use_dims)

set(T "proj_eb")
aphros_add_object(T ${T} proj_eb.cpp)
target_link_libraries(${T} PRIVATE utilconvdiff PUBLIC use_dims)

set(T "fluid_dummy")
aphros_add_object(T ${T} fluid_dummy.cpp)
target_link_libraries(${T} vars suspender posthook_default approx_eb)

set(T "electro")
aphros_add_object(T ${T} electro.cpp)
target_link_libraries(${T} PRIVATE approx_eb PUBLIC use_dims use_mpi)

