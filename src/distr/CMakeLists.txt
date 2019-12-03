set(T "local")
add_library(${T} local.cpp)
target_link_libraries(${T}
    hypre suspender vars parser report sysinfo openmp)

set(T "cubism")
add_library(${T} cubism.cpp)
target_link_libraries(${T} 
    hypre suspender vars parser report sysinfo hdf openmp)
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
target_compile_definitions(${T} PUBLIC _USE_HDF_)

set(T "cubismnc")
add_library(${T} cubismnc.cpp)
target_link_libraries(${T} 
    hypre suspender vars parser report sysinfo hdf openmp histogram)
target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
target_compile_definitions(${T} PUBLIC _USE_HDF_)

set(T "report")
add_library(${T} report.cpp)
target_link_libraries(${T})

set(T "distr")
add_library(${T} distrsolver.cpp)
target_link_libraries(${T} 
    parser suspender vars cubism cubismnc local dumper git subcomm hypresub)

set(T "distrbasic")
add_library(${T} distrbasic.cpp)
target_link_libraries(${T} distr)


