cmake_minimum_required(VERSION 3.12.0)

project("mfer")

#Require out - of - source build
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

option(USE_OPENMP "Enable OpenMP" ON)
option(USE_AVX "Enable AVX vectorized kernels. Requires a CPU with AVX2 and FMA support" ON)
option(USE_HYPRE "Enable Hypre for linear systems" ON)
option(USE_HDF "Enable HDF5 for output. Requres parallel HDF5 and MPI." ON)
option(USE_AMGX "Enable AMGX for linear systems, requires GPU with CUDA" OFF)
option(USE_MPI "Enable MPI for distributed computation" ON)
option(USE_FPZIP "Enable fpzip compression for communication (experimental)" OFF)
option(USE_OVERLAP "Enable Overlap library for exact initialization of ellipsoids. Requires Eigen3" OFF)
option(USE_DIM1 "Enable 1D geometry" OFF)
option(USE_DIM2 "Enable 2D geometry" OFF)
option(USE_DIM3 "Enable 3D geometry" ON)
option(USE_DIM4 "Enable 4D geometry" OFF)
option(USE_BACKEND_CUBISM "Enable Cubism backend for distributed computation. Requires MPI" ON)
option(USE_BACKEND_LOCAL "Enable local backend for shared-memory computation." ON)
option(USE_BACKEND_NATIVE "Enable native backend for distributed computation. Requires MPI" ON)

set(CMAKE_MACOSX_RPATH OFF)
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON)
set(CMAKE_LINK_DEPENDS_NO_SHARED ON)

set(APHROS_USE_MPI ${USE_MPI})
set(APHROS_USE_OPENMP ${USE_OPENMP})
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

if (USE_FPZIP)
  find_package(FPZIP REQUIRED)
endif()

if (USE_AMGX)
  find_package(CUDAToolkit)
endif()

if (USE_HYPRE)
  set(T "hypreext")
  set(HYPRE_DIR ${APHROS_PREFIX})
  add_library(${T} INTERFACE)
  target_include_directories(${T} INTERFACE ${HYPRE_DIR}/include)
  target_link_libraries(${T} INTERFACE -L${HYPRE_DIR}/lib -lHYPRE -lm)
  install(TARGETS ${T} EXPORT export)
endif()

if (USE_HDF)
  set(T "hdf")
  if (USE_MPI)
    set(HDF5_PREFER_PARALLEL ON)
  endif()
  set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
  find_package(HDF5 REQUIRED COMPONENTS C HL)
  if (USE_MPI AND NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "No parallel HDF5")
  endif()
  add_library(${T} INTERFACE)
  target_include_directories(${T} INTERFACE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${T} INTERFACE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  install(TARGETS ${T} EXPORT export)
endif()

set(T "use_avx")
add_library(${T} INTERFACE)
if (USE_AVX)
  target_compile_options(${T} INTERFACE -march=native)
endif ()
target_compile_definitions(${T} INTERFACE _USE_AVX_=$<BOOL:${USE_AVX}>)
install(TARGETS ${T} EXPORT export)

set(T "use_dims")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_DIM1_=$<BOOL:${USE_DIM1}>)
target_compile_definitions(${T} INTERFACE _USE_DIM2_=$<BOOL:${USE_DIM2}>)
target_compile_definitions(${T} INTERFACE _USE_DIM3_=$<BOOL:${USE_DIM3}>)
target_compile_definitions(${T} INTERFACE _USE_DIM4_=$<BOOL:${USE_DIM4}>)
install(TARGETS ${T} EXPORT export)

set(T "use_mpi")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_MPI_=$<BOOL:${USE_MPI}>)
install(TARGETS ${T} EXPORT export)

set(T "use_overlap")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_OVERLAP_=$<BOOL:${USE_OVERLAP}>)
install(TARGETS ${T} EXPORT export)

if (USE_HDF AND NOT USE_MPI)
  message(FATAL_ERROR "USE_HDF5 requires USE_MPI")
endif()
set(T "use_hdf")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_HDF_=$<BOOL:${USE_HDF}>)
install(TARGETS ${T} EXPORT export)

set(T "use_backend")
add_library(${T} INTERFACE)
target_compile_definitions(
  ${T} INTERFACE _USE_BACKEND_CUBISM_=$<BOOL:${USE_BACKEND_CUBISM}>)
target_compile_definitions(
  ${T} INTERFACE _USE_BACKEND_LOCAL_=$<BOOL:${USE_BACKEND_LOCAL}>)
target_compile_definitions(
  ${T} INTERFACE _USE_BACKEND_NATIVE_=$<BOOL:${USE_BACKEND_NATIVE}>)
install(TARGETS ${T} EXPORT export)

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
  convdiffe
  convdiffi
  convdiffvg
  curv
  distr
  distrbasic
  distrsolver
  dump
  dump_hdf
  dump_raw
  dump_xmf
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
  march
  mesh
  normal
  parser
  parse_template
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
  utilmpi
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

if (USE_BACKEND_CUBISM)
  set(ObjLibs ${ObjLibs} cubismnc)
endif()

if (USE_BACKEND_LOCAL)
  set(ObjLibs ${ObjLibs} local)
endif()

if (USE_BACKEND_NATIVE)
  set(ObjLibs ${ObjLibs} native comm_manager)
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
target_link_libraries(${EXE} aphros use_dims)

# test
add_subdirectory(test)
enable_testing()

configure_file(vars.cmake.in vars.cmake @ONLY)

install(EXPORT export FILE export.cmake DESTINATION cmake)
install(FILES client.cmake ${PROJECT_BINARY_DIR}/vars.cmake aphrosConfig.cmake DESTINATION cmake)
include(header.cmake)
