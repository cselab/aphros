cmake_minimum_required(VERSION 3.3.0)

project("mfer")

set(EXE "mfer")

set(CHPREFIX $ENV{CHPREFIX})

include("client.cmake")

# git revision
set(P "${PROJECT_SOURCE_DIR}/..")
set(C "${P}/src/util/gitgen.cpp")
# create on configure
if (NOT EXISTS ${C})
  execute_process(
    COMMAND "${P}/tool/gitrev" "${C}"
    WORKING_DIRECTORY "${P}"
  )
endif()
# create on build
add_custom_target(gitrev
  COMMAND "${P}/tool/gitrev" "${C}"
  WORKING_DIRECTORY "${P}"
)

include_directories(.)

function(add_library name)
  _add_library(${name} ${ARGN})
  install(TARGETS ${name} EXPORT export DESTINATION ${CHPREFIX}/lib)
endfunction(add_library)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(solver)
add_subdirectory(parse)
add_subdirectory(dump)
add_subdirectory(overlap)
add_subdirectory(young)

# hydro
set(T "hydro")
add_library(${T} INTERFACE)
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    "${CHPREFIX}/include/hydro")
set_property(TARGET ${T} APPEND PROPERTY INTERFACE_LINK_LIBRARIES 
    distr solver vof tvd convdiffvi simple overlap young utilfluid events)

# binary
add_executable(${EXE} main.cpp)
install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")
set_target_properties(${EXE} PROPERTIES PREFIX "ch.")
target_link_libraries(${EXE} hydro)

# test
add_subdirectory(test)
enable_testing()

install(EXPORT export FILE export.cmake DESTINATION ${CHPREFIX}/include)
install(FILES client.cmake hydro.cmake DESTINATION ${CHPREFIX}/include)
include(header.cmake)
