set(TT color column gyration filter split)

foreach(T ${TT})
  add_executable(${T} ${T}.c)
  target_link_libraries(${T} m csv table vtk chmath)
  install(
      TARGETS ${T}
      RENAME ${NAMEPREFIX}${T}
      DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endforeach()

