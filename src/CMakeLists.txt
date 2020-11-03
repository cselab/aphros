cmake_minimum_required(VERSION 3.12.0)

project("mfer")

#Require out - of - source build
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

option(USE_OPENMP "OpenMP support" ON)
option(USE_AVX "AVX vectorized kernels. Requires a CPU with AVX2 and FMA support" ON)
option(USE_HYPRE "Hypre for linear systems" ON)
option(USE_FPZIP "fpzip compression for communication (experimental)" OFF)

set(CMAKE_MACOSX_RPATH OFF)
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

include(client.cmake)
find_package(aphros/march)
find_package(aphros/inside)

# git revision
set(P "${PROJECT_SOURCE_DIR}/..")
set(C "${P}/src/util/gitgen.cpp")
# create on configure
if (NOT EXISTS ${C})
  execute_process(
    COMMAND "${P}/src/gitrev" "${C}"
    WORKING_DIRECTORY "${P}"
  )
endif()
# create on build
add_custom_target(gitrev
  COMMAND "${P}/src/gitrev" "${C}"
  WORKING_DIRECTORY "${P}"
)

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON)

include_directories(.)

# Override add_library().
# By default, creates OBJECT libraries.
function(add_library name)
  if ("${ARGN}" MATCHES "SHARED|STATIC|INTERFACE") # normal library
    _add_library(${name} ${ARGN})
    install(TARGETS ${name} EXPORT export DESTINATION lib)
  else() # object library
    _add_library(${name} OBJECT ${ARGN})
    install(TARGETS ${name} EXPORT export)
  endif()
  if (NOT "${ARGN}" MATCHES "INTERFACE")
    target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic -Wshadow)
  endif()
  target_include_directories(${name} INTERFACE "${APHROS_PREFIX}/include/aphros")
endfunction(add_library)

add_subdirectory(color)
add_subdirectory(distr)
add_subdirectory(dump)
add_subdirectory(func)
add_subdirectory(geom)
add_subdirectory(linear)
add_subdirectory(overlap)
add_subdirectory(parse)
add_subdirectory(solver)
add_subdirectory(util)
add_subdirectory(young)


set(ObjLibs
  allocator
  approx
  approx_eb
  argparse
  color
  convdiffe
  convdiffi
  convdiffvg
  cubismnc
  curv
  distr
  distrbasic
  distrsolver
  dump
  dump_hdf
  dumper
  electro
  embed
  events
  filesystem
  fluid_dummy
  format
  git
  gitgen
  histogram
  init
  init_contang
  init_vel
  linear
  local
  mesh
  normal
  overlap
  parser
  particles
  partstrmeshm
  posthook_default
  proj
  proj_eb
  report
  simple
  solver
  subcomm
  suspender
  sysinfo
  timer
  tracer
  utilconvdiff
  utildistr
  utilfluid
  utilhydro
  utillinear
  utilvof
  vars
  vof
  vofm
  young
)

if (USE_HYPRE)
  set(ObjLibs ${ObjLibs}
    hypre
    hypresub
    linear_hypre
  )
endif()

foreach(name ${ObjLibs})
  list(APPEND ObjLibsTarget "$<TARGET_OBJECTS:${name}>")
endforeach()

set(T "aphros")
add_library(${T} SHARED ${ObjLibsTarget})
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    "${APHROS_PREFIX}/include/aphros")
target_link_libraries(${T} ${ObjLibs})

# binary
set(EXE "mfer")
add_executable(${EXE} main.cpp)
install(TARGETS ${EXE} DESTINATION bin)
set_target_properties(${EXE} PROPERTIES PREFIX "ap.")
target_link_libraries(${EXE} aphros)

# test
add_subdirectory(test)
enable_testing()

install(EXPORT export FILE export.cmake DESTINATION cmake)
install(FILES client.cmake aphrosConfig.cmake DESTINATION cmake)
include(header.cmake)
