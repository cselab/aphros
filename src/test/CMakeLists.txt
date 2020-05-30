set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")

# find all directories containing CMakeLists.txt
set(l "")
file(GLOB g ${CSD}/*)
foreach(e ${g})
  if (EXISTS ${e}/CMakeLists.txt)
    LIST(APPEND l ${e})
  endif()
endforeach()

enable_testing()

# add them to project
foreach(e ${l})
  add_subdirectory(${e} ${e})
endforeach()

