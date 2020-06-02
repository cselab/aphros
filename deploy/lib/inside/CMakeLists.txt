set(T inside)
add_library(${T} SHARED
    bbox.c
    bbox.h
    err.c
    err.h
    inside.h
    main.c
    memory.c
    memory.h
    off.c
    ply.c
    predicate.c
    predicate.h
    predicate/main.c
    stl.c
  )
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "inside.h")
install(TARGETS ${T}
    EXPORT export DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
target_include_directories(${T} PUBLIC "${CMAKE_INSTALL_PREFIX}/include")
target_link_libraries(${T} PUBLIC m)
install(EXPORT export FILE insideConfig.cmake DESTINATION ${CMAKEPREFIX})
