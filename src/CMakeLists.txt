cmake_minimum_required(VERSION 3.12.0)

project("mfer")

#Require out - of - source build
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

option(USE_OPENMP "OpenMP support" ON)
option(USE_AVX "Enable AVX vectorized kernels. Requires a CPU with AVX2 and FMA support" ON)
option(USE_HYPRE "Enable Hypre for linear systems" ON)
option(USE_AMGX "Enable AMGX for linear systems, requires GPU with CUDA" OFF)
option(USE_FPZIP "Enable fpzip compression for communication (experimental)" OFF)
option(USE_OVERLAP "Enable Overlap library for exact initialization of ellipsoids. Requires Eigen3" ON)
option(USE_DIM1 "Enable 1D geometry" OFF)
option(USE_DIM2 "Enable 2D geometry" OFF)
option(USE_DIM3 "Enable 3D geometry" ON)
option(USE_DIM4 "Enable 4D geometry" OFF)

set(CMAKE_MACOSX_RPATH OFF)
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

include(client.cmake)

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

set(T "use_avx")
add_library(${T} INTERFACE)
if (USE_AVX)
  target_compile_definitions(${T} INTERFACE _USE_AVX_=1)
  target_compile_options(${T} INTERFACE -march=native)
else ()
  target_compile_definitions(${T} INTERFACE _USE_AVX_=0)
endif ()
install(TARGETS ${T} EXPORT export)

set(T "use_dims")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_DIM1_=$<BOOL:${USE_DIM1}>)
target_compile_definitions(${T} INTERFACE _USE_DIM2_=$<BOOL:${USE_DIM2}>)
target_compile_definitions(${T} INTERFACE _USE_DIM3_=$<BOOL:${USE_DIM3}>)
target_compile_definitions(${T} INTERFACE _USE_DIM4_=$<BOOL:${USE_DIM4}>)
install(TARGETS ${T} EXPORT export)

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
add_subdirectory(inside)
add_subdirectory(linear)
add_subdirectory(march)
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
  comm_manager
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
  dump_raw
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
  march
  mesh
  native
  normal
  parser
  particles
  partstrmeshm
  posthook_default
  primlist
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
set(BinLibs
  inside
)

if (USE_HYPRE)
  set(ObjLibs ${ObjLibs} hypre hypresub linear_hypre)
endif()

if (USE_AMGX)
  set(ObjLibs ${ObjLibs} linear_amgx)
endif()

if (USE_OVERLAP)
  set(ObjLibs ${ObjLibs} overlap)
endif()

foreach(name ${ObjLibs})
  list(APPEND ObjLibsTarget "$<TARGET_OBJECTS:${name}>")
endforeach()

set(T "aphros")
add_library(${T} SHARED ${ObjLibsTarget})
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    "${APHROS_PREFIX}/include/aphros")
target_link_libraries(${T} ${ObjLibs} ${BinLibs})

# binary
set(EXE "mfer")
add_executable(${EXE} main.cpp)
install(TARGETS ${EXE} DESTINATION bin)
set_target_properties(${EXE} PROPERTIES PREFIX "ap.")
target_link_libraries(${EXE} aphros)

# test
add_subdirectory(test)
enable_testing()

configure_file(vars.cmake.in vars.cmake @ONLY)

install(EXPORT export FILE export.cmake DESTINATION cmake)
install(FILES client.cmake ${PROJECT_BINARY_DIR}/vars.cmake aphrosConfig.cmake DESTINATION cmake)
include(header.cmake)
