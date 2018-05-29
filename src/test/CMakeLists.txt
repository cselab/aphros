set(r ${CMAKE_CURRENT_SOURCE_DIR})

# find all directories containing CMakeLists.txt
set(l "")
file(GLOB g ${r}/*)
foreach(e ${g})
  if (EXISTS ${e}/CMakeLists.txt)
    LIST(APPEND l ${e})
  endif()
endforeach()

# add them to project
foreach(e ${l})
  add_subdirectory(${e})
endforeach()

enable_testing()
