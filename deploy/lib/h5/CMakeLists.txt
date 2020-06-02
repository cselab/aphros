include(h5.dep.cmake)

set(T h5)
set(H h5.h h5serial.h)
add_library(${T} SHARED main.c h5.h h5serial.h)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "${H}")
target_link_libraries(${T} PUBLIC hdf)
target_include_directories(${T} PUBLIC "${CMAKE_INSTALL_PREFIX}/include")
install(TARGETS ${T}
    EXPORT export_${T} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(EXPORT export_${T} FILE h5.lib.cmake DESTINATION ${CMAKEPREFIX})
install(FILES h5.dep.cmake h5Config.cmake DESTINATION ${CMAKEPREFIX})
