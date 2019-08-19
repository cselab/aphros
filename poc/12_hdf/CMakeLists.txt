cmake_minimum_required(VERSION 3.3.0)
set(CHPREFIX $ENV{CHPREFIX})
project("hdf" C)

include(h5.dep.cmake)

set(T h5)
add_library(${T} SHARED h5.c h5.h)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER h5.h)
target_link_libraries(${T} PUBLIC hdf)
install(TARGETS ${T} 
    EXPORT export DESTINATION ${CHPREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CHPREFIX}/include)
install(EXPORT export FILE h5.lib.cmake DESTINATION ${CHPREFIX}/include)
install(FILES h5.dep.cmake h5.cmake DESTINATION ${CHPREFIX}/include)
