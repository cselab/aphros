set(TARGL "local")
add_library(${TARGL} local.cpp)
target_link_libraries(${TARGL} hypre suspender vars interp)


set(TARGC "cubism")
add_library(${TARGC} cubism.cpp)
target_link_libraries(${TARGC} hypre suspender vars interp)
# cubism specific
target_compile_definitions(${TARGC} PUBLIC _ALIGNBYTES_=16)
# HDF5
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
target_include_directories(${TARGC} PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(${TARGC} ${HDF5_LIBRARIES})
target_compile_definitions(${TARGC} PUBLIC _USE_HDF_)


set(TARG "distr")
add_library(${TARG} distrsolver.cpp)
target_link_libraries(${TARG} interp suspender vars cubism local)

