cmake_minimum_required(VERSION 3.10.2)

project("mfer")

# Allow using <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

# Require out - of - source build
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

option(USE_AMGX "Enable AMGX for linear systems, requires GPU with CUDA" OFF)
option(USE_AVX "Enable AVX vectorized kernels. Requires a CPU with AVX2 and FMA support" ON)
option(USE_BACKEND_CUBISM "Enable Cubism backend for distributed computation. Requires MPI" ON)
option(USE_BACKEND_LOCAL "Enable local backend for shared-memory computation." ON)
option(USE_BACKEND_NATIVE "Enable native backend for distributed computation. Requires MPI" ON)
option(USE_DIM1 "Enable 1D geometry" OFF)
option(USE_DIM2 "Enable 2D geometry" ON)
option(USE_DIM3 "Enable 3D geometry" ON)
option(USE_DIM4 "Enable 4D geometry" OFF)
option(USE_FPZIP "Enable fpzip compression for communication (experimental)" OFF)
option(USE_HDF "Enable HDF5 for output. Requres parallel HDF5 and MPI." ON)
option(USE_HYPRE "Enable Hypre for linear systems" ON)
option(USE_MPI "Enable MPI for distributed computation" ON)
option(USE_OPENMP "Enable OpenMP" ON)
option(USE_OPENCL "Enable OpenCL for GPU support" OFF)
option(USE_TESTS "Build tests and enable testing" ON)
option(USE_WARNINGS "Enable compiler warnings" ON)

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
  install(TARGETS ${T} EXPORT export DESTINATION lib)
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

set(T "use_opencl")
add_library(${T} INTERFACE)
target_compile_definitions(${T} INTERFACE _USE_OPENCL_=$<BOOL:${USE_OPENCL}>)
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

# Defines an object library named ${T}_obj and an interface library named ${T}.
function(add_object T)
  add_library(${T}_obj OBJECT ${ARGN})
  add_library(${T} INTERFACE)
  install(TARGETS ${T} EXPORT export DESTINATION lib)
  #install(TARGETS ${T}_obj EXPORT export DESTINATION lib)
  if (USE_WARNINGS)
    target_compile_options(${T}_obj PRIVATE -Wall -Wextra -pedantic -Wshadow)
  endif()
endfunction(add_object)

function(object_link_libraries T)
  cmake_parse_arguments(L "" "" "PRIVATE;PUBLIC;INTERFACE" ${ARGN})

  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_UNPARSED_ARGUMENTS})
    target_link_libraries(${T} INTERFACE ${L})
  endforeach()

  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_PRIVATE} ${L_UNPARSED_ARGUMENTS})
    target_compile_definitions(${T}_obj
      PRIVATE $<TARGET_PROPERTY:${L},INTERFACE_COMPILE_DEFINITIONS>)
    target_compile_options(${T}_obj
      PRIVATE $<TARGET_PROPERTY:${L},INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(${T}_obj
      PRIVATE $<TARGET_PROPERTY:${L},INTERFACE_INCLUDE_DIRECTORIES>)
  endforeach()
endfunction()

function(object_compile_definitions T)
  cmake_parse_arguments(L "" "" "PRIVATE;PUBLIC;INTERFACE" ${ARGN})

  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_UNPARSED_ARGUMENTS})
    target_compile_definitions(${T} INTERFACE ${L})
  endforeach()

  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_PRIVATE} ${L_UNPARSED_ARGUMENTS})
    target_compile_definitions(${T}_obj PRIVATE ${L})
  endforeach()
endfunction()

function(object_compile_options T)
  cmake_parse_arguments(L "" "" "PRIVATE;PUBLIC;INTERFACE" ${ARGN})
  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_UNPARSED_ARGUMENTS})
    target_compile_options(${T} INTERFACE ${L})
  endforeach()
  foreach (L ${L_INTERFACE} ${L_PUBLIC} ${L_PRIVATE} ${L_UNPARSED_ARGUMENTS})
    target_compile_options(${T}_obj PRIVATE ${L})
  endforeach()
endfunction()

function(get_link_libraries0 T)
  list(APPEND VISITED ${T})
  get_target_property(LINKED ${T} INTERFACE_LINK_LIBRARIES)
  foreach(L ${LINKED})
    if (TARGET ${L})
      list(FIND VISITED ${L} V)
      if (${V} EQUAL -1)
        get_link_libraries0(${L})
        list(APPEND RES ${L})
      endif()
    endif()
  endforeach()
  set(VISITED "${VISITED}" PARENT_SCOPE)
  set(RES "${RES}" PARENT_SCOPE)
endfunction()

function(get_link_libraries VARNAME T)
  set(VISITED)
  set(RES)
  get_link_libraries0(${T})
  set(${VARNAME} "${RES}" PARENT_SCOPE)
endfunction()

function(add_executable_with_objects T src)
  set(LIBS ${ARGN})
  get_object_files(P ${LIBS})
  add_executable(${T} ${src} ${P})
  target_link_libraries(${T} ${LIBS})
endfunction()

function(get_object_files VARNAME)
  set(RES)
  foreach (SRC ${ARGN})
    get_link_libraries(LINKED ${SRC})
    list(APPEND LINKED "${SRC}")
    foreach (L ${LINKED})
      if (TARGET ${L}_obj)
        list(APPEND RES "$<TARGET_OBJECTS:${L}_obj>")
      endif()
    endforeach()
  endforeach()
  set(${VARNAME} ${RES} PARENT_SCOPE)
endfunction()


add_subdirectory(color)
add_subdirectory(distr)
add_subdirectory(dump)
add_subdirectory(func)
add_subdirectory(geom)
add_subdirectory(inside)
add_subdirectory(linear)
add_subdirectory(march)
add_subdirectory(opencl)
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
  hydro
  hydro_post
  init
  init_contang
  init_vel
  linear
  logger
  march
  mesh
  normal
  overlap
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
  utilsystem
  utilvof
  vars
  visual
  vof
  vofm
  young
)
set(UseLibs
  use_dims
  use_mpi
  use_backend
  use_overlap
)
set(StaticLibs
  inside
)

if (USE_HYPRE)
  set(ObjLibs ${ObjLibs} hypre hypresub linear_hypre)
endif()

if (USE_AMGX)
  set(ObjLibs ${ObjLibs} linear_amgx)
endif()

if (USE_OPENCL)
  set(ObjLibs ${ObjLibs} opencl linear_conjugate_cl)
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

set(ObjLibsTarget)
foreach(name ${ObjLibs})
  list(APPEND ObjLibsTarget "$<TARGET_OBJECTS:${name}_obj>")
endforeach()

set(T aphros)
add_library(${T} SHARED ${ObjLibsTarget})
target_include_directories(${T} INTERFACE "${APHROS_PREFIX}/include/aphros")
target_link_libraries(${T} PUBLIC ${ObjLibs} ${StaticLibs} ${UseLibs} openmp)
install(TARGETS ${T} EXPORT export DESTINATION lib)

add_subdirectory(kernel)

set(T aphros_c)
add_library(${T} SHARED aphros_c/git.cpp aphros_c/main.cpp aphros_c/parser.cpp)
target_link_libraries(${T} PUBLIC aphros ${UseLibs})
install(TARGETS ${T} EXPORT export DESTINATION lib)

# executables

set(T conf2py)
add_executable_with_objects(${T} parse/conf2py.cpp parser format argparse)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
install(TARGETS ${T} EXPORT export DESTINATION bin)

set(T mfer)
add_executable(${T} main.c)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
target_link_libraries(${T} aphros_c)
install(TARGETS ${T} EXPORT export DESTINATION bin)

add_custom_command(
  OUTPUT  explorer.inc
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/explorer_gen"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_custom_target(explorer_inc DEPENDS explorer.inc)

set(T explorer)
add_executable(${T} explorer.cpp)
add_dependencies(${T} explorer_inc)
target_link_libraries(${T} aphros)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
install(TARGETS ${T} EXPORT export DESTINATION bin)


if (USE_TESTS)
  add_subdirectory(test)
  enable_testing()
endif()

configure_file(vars.cmake.in vars.cmake @ONLY)

install(EXPORT export FILE export.cmake DESTINATION cmake)
install(FILES client.cmake ${PROJECT_BINARY_DIR}/vars.cmake aphrosConfig.cmake DESTINATION cmake)
include(header.cmake)
