set(T vof)
add_library(${T} SHARED vof.h vof.c)
target_link_libraries(${T} vofi)
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "vof.h")
install(TARGETS ${T}
    EXPORT export_${T} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(EXPORT export_${T} FILE vofConfig.cmake DESTINATION ${CMAKEPREFIX})
