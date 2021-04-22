set(T suspender)
add_object(${T} suspender.cpp)

set(T "utilsystem")
add_object(${T} system.c)

set(T "sysinfo")
add_object(${T} sysinfo.cpp)
object_link_libraries(${T} utilmpi utilsystem PRIVATE openmp)
object_compile_definitions(${T} PRIVATE _USE_AMGX_=$<BOOL:${USE_AMGX}>)

set(T "utilhydro")
add_object(${T} hydro.cpp)
object_link_libraries(${T} vars dumper approx vof vofm use_dims)

set(T "utillinear")
add_object(${T} linear.cpp)
object_link_libraries(${T} vars linear)

set(T "utilfluid")
add_object(${T} fluid.cpp)
object_link_libraries(${T} approx_eb)

set(T "visual")
add_object(${T} visual.cpp)
object_link_libraries(${T} mesh)

set(T "utilconvdiff")
add_object(${T} convdiff.cpp)
object_link_libraries(${T} convdiffvg)

set(T "utilvof")
add_object(${T} vof.cpp)
object_link_libraries(${T} march approx)

set(T "utildistr")
add_object(${T} distr.cpp)
object_link_libraries(${T} use_mpi)

set(T "utilmpi")
add_object(${T} mpi.cpp)
object_link_libraries(${T} use_mpi)

set(T "events")
add_object(${T} events.cpp)
object_link_libraries(${T} vars parser)

set(T "timer")
add_object(${T} timer.cpp)

set(T gitgen)
add_object(${T} gitgen.cpp)
add_dependencies(${T} gitrev)

set(T "filesystem")
add_object(${T} filesystem.cpp)
object_link_libraries(${T} utilsystem)

set(T "git")
add_object(${T} git.cpp)
object_link_libraries(${T} gitgen)

set(T "posthook_default")
add_object(${T} posthook_default.cpp)
object_link_libraries(${T} suspender vars use_mpi use_dims)

set(T "subcomm")
if (USE_OPENMP AND USE_MPI)
  add_object(${T} subcomm.cpp)
else()
  add_object(${T} subcomm_dummy.cpp)
endif()
object_link_libraries(${T} sysinfo utilsystem PRIVATE openmp)

set(T "histogram")
add_object(${T} histogram.cpp)
object_link_libraries(${T} timer logger use_mpi)

set(T "logger")
add_object(${T} logger.cpp)

set(T "format")
add_object(${T} format.cpp)
object_link_libraries(${T} logger)

set(T "allocator")
add_object(${T} fixed_allocator.cpp)

set(T hydro_post)
add_object(${T} hydro_post.cpp)
object_link_libraries(${T} PUBLIC dump visual utilhydro PRIVATE use_dims use_mpi)
