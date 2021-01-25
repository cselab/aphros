set(T csv)
add_library(${T} SHARED
    main.c
    memory.c
    err.c
  )
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "${T}.h")
install(TARGETS ${T} 
    EXPORT export_${T} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
target_include_directories(${T} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(${T} PUBLIC m)
install(EXPORT export_${T} FILE ${T}Config.cmake DESTINATION ${CMAKEPREFIX})
