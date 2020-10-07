cmake_minimum_required(VERSION 3.3.0)

project("mfer")

#Require out - of - source build
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
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

set(ObjLibs)

# Override add_library().
# By default, creates OBJECT libraries.
function(add_library name)
  if ("${ARGN}" MATCHES "SHARED|STATIC|INTERFACE")
    _add_library(${name} ${ARGN})
  else()
    set(ObjLibs "${ObjLibs}" "${name}" CACHE INTERNAL "ObjLibs")
    _add_library(${name} OBJECT ${ARGN})
  endif()
  if (NOT "${ARGN}" MATCHES "INTERFACE")
    target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic -Wshadow)
  endif()
  install(TARGETS ${name} EXPORT export DESTINATION lib)
  target_include_directories(${name} INTERFACE "${APHROS_PREFIX}/include/aphros")
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


set(ObjLibsTarget)
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
