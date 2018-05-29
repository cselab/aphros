set(TARG "distr")

include_directories(.)
add_library(${TARG} 
  distrsolver.cpp 
  local.cpp 
  cubism.cpp 
  ../util/suspender.cpp 
  ../dump/dumper.cpp
  ../parse/vars.cpp
  ../parse/interp.cpp
  )

# OpenMP support                              
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              


# HDF5
set(HDF5_PREFER_PARALLEL on)
find_package(HDF5 REQUIRED)
add_definitions(" -D_USE_HDF_")
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(${TARG} ${HDF5_LIBRARIES})

# cubism specific
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALIGNBYTES_=16") 

add_subdirectory(../linear)
target_link_libraries(${TARG} hypre)

