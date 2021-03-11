add_library(suspender suspender.cpp)

set(T "sysinfo")
add_library(${T} sysinfo.cpp)
target_link_libraries(${T} openmp utilmpi)
target_compile_definitions(${T} PRIVATE _USE_AMGX_=$<BOOL:${USE_AMGX}>)

set(T "utilhydro")
add_library(${T} hydro.cpp)
target_link_libraries(${T} PUBLIC vars dumper approx vof vofm PRIVATE use_dims)

set(T "utillinear")
add_library(${T} linear.cpp)
target_link_libraries(${T} vars linear)

set(T "utilfluid")
add_library(${T} fluid.cpp)
target_link_libraries(${T})

set(T "utilconvdiff")
add_library(${T} convdiff.cpp)
target_link_libraries(${T} convdiffvg)

set(T "utilvof")
add_library(${T} vof.cpp)
target_link_libraries(${T} march approx)

set(T "utildistr")
add_library(${T} distr.cpp)
target_link_libraries(${T} use_mpi)

set(T "utilmpi")
add_library(${T} mpi.cpp)
target_link_libraries(${T} use_mpi)

set(T "events")
add_library(${T} events.cpp)
target_link_libraries(${T} vars parser)

set(T "timer")
add_library(${T} timer.cpp)

add_library(gitgen gitgen.cpp)
add_dependencies(gitgen gitrev)

set(T "filesystem")
add_library(${T} filesystem.cpp)

set(T "git")
add_library(${T} git.cpp)
target_link_libraries(${T} gitgen)

set(T "posthook_default")
add_library(${T} posthook_default.cpp)
target_link_libraries(${T} suspender vars use_mpi use_dims)

set(T "subcomm")
if (USE_OPENMP AND USE_MPI)
  add_library(${T} subcomm.cpp)
else()
  add_library(${T} subcomm_dummy.cpp)
endif()
target_link_libraries(${T} sysinfo openmp)

set(T "histogram")
add_library(${T} histogram.cpp)
target_link_libraries(${T} timer logger use_mpi)

set(T "logger")
add_library(${T} logger.cpp)

set(T "format")
add_library(${T} format.cpp)
target_link_libraries(${T} logger)

set(T "system")
add_library(${T} system.cpp)

set(T "allocator")
add_library(${T} fixed_allocator.cpp)
