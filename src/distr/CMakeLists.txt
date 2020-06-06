set(T "local")
add_library(${T} local.cpp)
target_link_libraries(${T} distr)

set(T "cubism")
add_library(${T} cubism.cpp)
target_link_libraries(${T} distr hdf)
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
target_compile_definitions(${T} PUBLIC _USE_HDF_)

set(T "cubismnc")
add_library(${T} cubismnc.cpp)
target_link_libraries(${T} hdf distr)
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
target_compile_definitions(${T} PUBLIC _USE_HDF_)
if (FPZIP_FOUND)
    target_link_libraries(${T} fpzip)
    target_compile_definitions(${T} PUBLIC _USE_FPZIP_)
endif ()

set(T "report")
add_library(${T} report.cpp)
target_link_libraries(${T})

set(T "distr")
add_library(${T} distr.cpp)
target_link_libraries(${T}
    sysinfo report hypre parser suspender vars histogram
    dumper git subcomm hypresub openmp)

set(T "distrsolver")
add_library(${T} distrsolver.cpp)
target_link_libraries(${T} distr local cubism cubismnc)

set(T "distrbasic")
add_library(${T} distrbasic.cpp)
target_link_libraries(${T} distrsolver)


