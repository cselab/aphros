set(T "local")
add_library(${T} local.cpp)
target_link_libraries(${T} distr)

set(T "cubismnc")
add_library(${T} cubismnc.cpp)
target_link_libraries(${T} hdf distr)
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
target_compile_definitions(${T} PUBLIC _USE_HDF_)
target_compile_options(${T} PRIVATE -Wall -Wextra -pedantic -Wno-shadow)
if (USE_FPZIP)
  target_link_libraries(${T} fpzip)
endif ()
target_compile_definitions(${T} PRIVATE _USE_FPZIP_=$<BOOL:${USE_FPZIP}>)

set(T "native")
add_library(${T} native.cpp)
target_link_libraries(${T} hdf distr)

set(T "report")
add_library(${T} report.cpp)
target_link_libraries(${T})

set(T "distr")
add_library(${T} distr.cpp)
target_link_libraries(${T}
    sysinfo report parser suspender vars histogram dumper git subcomm openmp)

set(T "distrsolver")
add_library(${T} distrsolver.cpp)
target_link_libraries(${T} distr local cubismnc native utildistr)
target_compile_definitions(
  ${T} PRIVATE _USE_BACKEND_CUBISM_=$<BOOL:${USE_BACKEND_CUBISM}>)
target_compile_definitions(
  ${T} PRIVATE _USE_BACKEND_LOCAL_=$<BOOL:${USE_BACKEND_LOCAL}>)
target_compile_definitions(
  ${T} PRIVATE _USE_BACKEND_NATIVE_=$<BOOL:${USE_BACKEND_NATIVE}>)

set(T "distrbasic")
add_library(${T} distrbasic.cpp)
target_link_libraries(${T} distrsolver utildistr)

set(T "comm_manager")
add_library(${T} comm_manager.cpp)
target_link_libraries(${T} distr)
