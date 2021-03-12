set(T suspender)
aphros_add_object(T ${T} suspender.cpp)

set(T "utilsystem")
aphros_add_object(T ${T} system.c)

set(T "sysinfo")
aphros_add_object(T ${T} sysinfo.cpp)
target_link_libraries(${T} openmp utilmpi utilsystem)
target_link_libraries(sysinfo openmp)
target_compile_definitions(${T} PRIVATE _USE_AMGX_=$<BOOL:${USE_AMGX}>)

set(T "utilhydro")
aphros_add_object(T ${T} hydro.cpp)
target_link_libraries(${T} PUBLIC vars dumper approx vof vofm PRIVATE use_dims)

set(T "utillinear")
aphros_add_object(T ${T} linear.cpp)
target_link_libraries(${T} vars linear)

set(T "utilfluid")
aphros_add_object(T ${T} fluid.cpp)
target_link_libraries(${T})

set(T "utilconvdiff")
aphros_add_object(T ${T} convdiff.cpp)
target_link_libraries(${T} convdiffvg)

set(T "utilvof")
aphros_add_object(T ${T} vof.cpp)
target_link_libraries(${T} march approx)

set(T "utildistr")
aphros_add_object(T ${T} distr.cpp)
target_link_libraries(${T} use_mpi)

set(T "utilmpi")
aphros_add_object(T ${T} mpi.cpp)
target_link_libraries(${T} use_mpi)

set(T "events")
aphros_add_object(T ${T} events.cpp)
target_link_libraries(${T} vars parser)

set(T "timer")
aphros_add_object(T ${T} timer.cpp)

set(T gitgen)
aphros_add_object(T ${T} gitgen.cpp)
add_dependencies(${T} gitrev)

set(T "filesystem")
aphros_add_object(T ${T} filesystem.cpp)
target_link_libraries(${T} utilsystem)

set(T "git")
aphros_add_object(T ${T} git.cpp)
target_link_libraries(${T} gitgen)

set(T "posthook_default")
aphros_add_object(T ${T} posthook_default.cpp)
target_link_libraries(${T} suspender vars use_mpi use_dims)

set(T "subcomm")
if (USE_OPENMP AND USE_MPI)
  aphros_add_object(T ${T} subcomm.cpp)
else()
  aphros_add_object(T ${T} subcomm_dummy.cpp)
endif()
target_link_libraries(${T} sysinfo openmp utilsystem)

set(T "histogram")
aphros_add_object(T ${T} histogram.cpp)
target_link_libraries(${T} timer logger use_mpi)

set(T "logger")
aphros_add_object(T ${T} logger.cpp)

set(T "format")
aphros_add_object(T ${T} format.cpp)
target_link_libraries(${T} logger)

set(T "allocator")
aphros_add_object(T ${T} fixed_allocator.cpp)
