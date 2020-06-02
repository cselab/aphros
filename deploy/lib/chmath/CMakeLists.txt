set(T chmath)
add_library(${T} SHARED
    main.c
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
find_package(GSL REQUIRED)
target_link_libraries(${T} PUBLIC m GSL::gsl)
install(EXPORT export_${T} FILE ${T}Config.cmake DESTINATION ${CMAKEPREFIX})
