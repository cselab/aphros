set(TARG "distr")

add_library(${TARG} 
  distrsolver.cpp 
  local.cpp 
  cubism.cpp 
  ../util/suspender.cpp 
  ../dump/dumper.cpp
  ../parse/vars.cpp
  ../parse/interp.cpp
  )

target_include_directories(${TARG} PUBLIC ../)

# HDF5
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
target_include_directories(${TARG} PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(${TARG} ${HDF5_LIBRARIES})
target_compile_definitions(${TARG} PUBLIC _USE_HDF_)

# cubism specific
target_compile_definitions(${TARG} PUBLIC _ALIGNBYTES_=16)

target_link_libraries(${TARG} hypre)


