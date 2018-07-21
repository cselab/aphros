set(T "local")
add_library(${T} local.cpp)
target_link_libraries(${T} hypre suspender vars interp report sysinfo)

set(T "cubism")
add_library(${T} cubism.cpp)
target_link_libraries(${T} hypre suspender vars interp report sysinfo)
# cubism specific
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
# HDF5
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
target_include_directories(${T} PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(${T} ${HDF5_LIBRARIES})
target_compile_definitions(${T} PUBLIC _USE_HDF_)

set(T "report")
add_library(${T} report.cpp)
target_link_libraries(${T})

set(T "distr")
add_library(${T} distrsolver.cpp)
target_link_libraries(${T} interp suspender vars cubism local dumper)

