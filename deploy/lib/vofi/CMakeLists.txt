set(T vofi)
add_library(${T} SHARED
    checkconsistency.c
    getcc.c
    getdirs.c
    getfh.c
    getintersections.c
    getlimits.c
    getmin.c
    getzero.c
    integrate.c
    interface.c
    vofi_GL.h
    vofi.h
    vofi_stddecl.h)
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "vofi.h")
install(TARGETS ${T}
    EXPORT export_${T} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
target_include_directories(${T} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)
include_directories(.)
target_link_libraries(${T} PUBLIC m)
install(EXPORT export_${T} FILE vofiConfig.cmake DESTINATION ${CMAKEPREFIX})
