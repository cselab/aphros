set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_LINK_DEPENDS_NO_SHARED TRUE)

# find all directories containing CMakeLists.txt
set(list "")
file(GLOB glob ${CSD}/*)
foreach(dir ${glob})
  if (EXISTS ${dir}/CMakeLists.txt)
    LIST(APPEND list ${dir})
  endif()
endforeach()

enable_testing()

# Adds test using the current source directory name as prefix
# NAME: test name suffix
# COMMAND: command to run
function(add_test_current)
  cmake_parse_arguments(CURPREFIX "" "NAME" "COMMAND" ${ARGN})
  get_filename_component(prefix ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  if (DEFINED CURPREFIX_NAME)
    set(name ${prefix}-${CURPREFIX_NAME})
  else()
    set(name ${prefix})
  endif()
  add_test(NAME ${name} COMMAND ${CURPREFIX_COMMAND})
  set_tests_properties(${name} PROPERTIES RESOURCE_LOCK ${prefix})
endfunction()


foreach(dir ${list})
  add_subdirectory(${dir} ${dir})
endforeach()

