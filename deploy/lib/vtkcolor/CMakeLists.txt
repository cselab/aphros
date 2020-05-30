set(T vtkcolor)
add_executable(${T} color.c)
target_link_libraries(${T} m)
install(
    TARGETS ${T}
    RENAME ${NAMEPREFIX}${T}
    DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
