cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

# MPI flag
if (DEFINED MPI)
  set(MPI_ENABLE ${MPI})
else()
  set(MPI_ENABLE ON)
endif()

# Static build flag
if (DEFINED STATIC)
  set(STATIC_BUILD ${STATIC})   
else()
  set(STATIC_BUILD OFF)
endif()

# Debug flag
#if (NOT DEFINED DEBUG)
#  set(DEBUG OFF)   
#endif()

#if (DEBUG)
#  set(CMAKE_BUILD_TYPE Debug)
#else()
#  set(CMAKE_BUILD_TYPE Release)
#endif()

project("mfer")

set(BINARY_NAME "mfer")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  

# OpenMP support                              
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS 
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()                              


# GCC warnings and optimizations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic ")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

if (UNIX)
  #if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic --std=c++11")
  #endif()
  #if (DEBUG)
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
  #else()
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")
  #  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")  # remove symbolic information
  #endif()
endif()

# Source folders
aux_source_directory(Cubism/ SOURCES_CUBISM)
aux_source_directory(CubismDistr/ SOURCES_CUBISM_DISTR)
aux_source_directory(hydro/ SOURCES_HYDRO)

# Include dirs
include_directories(SOUCRES_CUBISM)
include_directories(.)

set(SOURCE_FILES
  ${SOURCES_CUBISM}
  ${SOURCES_CUBISM_DISTR}
  ${SOURCES_HYDRO}
)

# Add target
add_executable(${BINARY_NAME} ${SOURCE_FILES})
     
# C++11
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BINARY_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# Static build
if (STATIC_BUILD)
  SET(BUILD_SHARED_LIBRARIES OFF)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_START_STATIC 1)
  set_property(TARGET ${BINARY_NAME} PROPERTY LINK_SEARCH_END_STATIC 1)
  if (UNIX)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  else()
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
  endif()
endif()

# MPI
if (MPI_ENABLE) 
  find_package(MPI REQUIRED)
  
  add_definitions("-DMPI_ENABLE")
  
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  
  #include_directories(${MPI_INCLUDE_PATH})
  #target_link_libraries(${BINARY_NAME} ${MPI_LIBRARIES})

  #if (MPI_COMPILE_FLAGS)
  #  set_target_properties(${BINARY_NAME} PROPERTIES
  #    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  #endif()

  #if (MPI_LINK_FLAGS)
  #  set_target_properties(${BINARY_NAME} PROPERTIES
  #    LINK_FLAGS "${MPI_LINK_FLAGS}")
  #endif()
endif()

# HDF5
if (1) 
  set(HDF5_PREFER_PARALLEL on)
  find_package(HDF5 REQUIRED)
  
  add_definitions(" -D_USE_HDF_")
  
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(${BINARY_NAME} ${HDF5_LIBRARIES})
endif()

# hypre
include_directories(${HYPRE_DIR}/include)
#link_directories(${HYPRE_DIR}/lib) # no effect
target_link_libraries(${BINARY_NAME} -lHYPRE -lm)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HYPRE_DIR}/lib")

# cubism specific
set(BS 8)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_BLOCKSIZE_=${BS} -D_BLOCKSIZEX_=${BS} -D_BLOCKSIZEY_=${BS} -D_BLOCKSIZEZ_=${BS}") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ALIGNBYTES_=16") 

target_link_libraries(${BINARY_NAME} ${LIBRARIES})

install(TARGETS ${BINARY_NAME} RUNTIME DESTINATION $ENV{HOME}/bin)


# tests
enable_testing()

add_test(t1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${BINARY_NAME} a b)
add_test(t2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${BINARY_NAME} a a)

