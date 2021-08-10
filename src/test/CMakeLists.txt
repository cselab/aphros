set(CMAKE_LINK_DEPENDS_NO_SHARED TRUE)

# Find all directories containing CMakeLists.txt
set(Tests "")
file(GLOB Glob ${CMAKE_CURRENT_SOURCE_DIR}/*)
foreach(Test ${Glob})
  if (EXISTS ${Test}/CMakeLists.txt)
    list(APPEND Tests ${Test})
  endif()
endforeach()

enable_testing()

# Adds test using the current source directory name as prefix
# NAME: test name suffix
# COMMAND: command to run
# ASSETS: list of asset files to copy from source to binary directory
#         defaults to all files in source directory
function(add_test_current)
  cmake_parse_arguments(RESULT "" "NAME" "COMMAND;ASSETS" ${ARGN})
  get_filename_component(prefix ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  if (DEFINED RESULT_NAME)
    set(name ${prefix}-${RESULT_NAME})
  else()
    set(name ${prefix})
  endif()
  add_test(NAME ${name} COMMAND ${RESULT_COMMAND})
  set_tests_properties(${name} PROPERTIES RESOURCE_LOCK ${prefix})
  copy_test_assets(${RESULT_ASSETS})
endfunction()


# Copies files specified in arguments to binary directory of current test
# If no arguments given, copies all files in current source directory
function(copy_test_assets)
  if (${ARGC} EQUAL 0)
    file(GLOB Assets *)
    list(FILTER Assets EXCLUDE REGEX ".*.gitignore$")
    list(FILTER Assets EXCLUDE REGEX ".*.cpp$")
    list(FILTER Assets EXCLUDE REGEX ".*.h$")
    list(FILTER Assets EXCLUDE REGEX ".*CMakeLists.txt$")
  else()
    set(Assets ${ARGN})
  endif()
  file(COPY ${Assets} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

foreach(Test ${Tests})
  add_subdirectory(${Test})
endforeach()


