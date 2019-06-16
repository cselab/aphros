add_library(suspender suspender.cpp)
add_library(sysinfo sysinfo.cpp)

set(T "utilfluid")
add_library(${T} fluid.cpp)
target_link_libraries(${T} vars)

set(T "events")
add_library(${T} events.cpp)
target_link_libraries(${T} vars parser)

set(T "timer")
add_library(${T} timer.cpp)

set(T "git")
add_library(${T} git.cpp)
target_link_libraries(${T} gitgen)

set(C "${CMAKE_CURRENT_SOURCE_DIR}/gitgen.cpp")
add_library(gitgen ${C})
add_custom_target(gitrev ALL DEPENDS ${C}_)
add_custom_command(
  OUTPUT ${C}_
  PRE_BUILD
  COMMAND "${PROJECT_SOURCE_DIR}/../tool/gitrev" ${C}
  COMMENT "Update Git revision"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
