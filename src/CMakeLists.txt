cmake_minimum_required(VERSION 3.3.0)

project("mfer")

# Disallow in-source build
if (CMAKE_SOURCE_DIR} STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

option(USE_OPENMP "Build with OpenMP support." ON)

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

function(add_library name)
  if ("${ARGN}" MATCHES "INTERFACE")
    _add_library(${name} ${ARGN})
  else()
    _add_library(${name} SHARED ${ARGN})
    set_target_properties(${name} PROPERTIES PREFIX "libaphros_")
    target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic -Wshadow)
  endif()
  install(TARGETS ${name} EXPORT export DESTINATION lib)
  target_include_directories(${name} INTERFACE "${APHROS_PREFIX}/include/hydro")
endfunction(add_library)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(func)
add_subdirectory(solver)
add_subdirectory(parse)
add_subdirectory(dump)
add_subdirectory(overlap)
add_subdirectory(young)
add_subdirectory(color)

# hydro
set(T "hydro")
add_library(${T} INTERFACE)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    "${APHROS_PREFIX}/include/hydro")
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES
    distrsolver solver vof vofm normal convdiffvg simple proj proj_eb
    fluid_dummy init init_vel init_contang
    utilhydro events curv approx posthook_default tracer
    particles dump_hdf electro linear utillinear distrbasic)

# binary
set(EXE "mfer")
add_executable(${EXE} main.cpp)
install(TARGETS ${EXE} DESTINATION bin)
set_target_properties(${EXE} PROPERTIES PREFIX "ap.")
target_link_libraries(${EXE} hydro)

# test
add_subdirectory(test)
enable_testing()

install(EXPORT export FILE export.cmake DESTINATION cmake/aphros)
install(FILES client.cmake hydroConfig.cmake DESTINATION cmake/aphros)
include(header.cmake)
