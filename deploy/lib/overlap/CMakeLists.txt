set(T overlap)
add_library(${T} SHARED
  overlap.cpp
  overlap.h
  overlap.inc
)
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "overlap.h")
install(TARGETS ${T}
    EXPORT export_${T} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
target_include_directories(${T} PUBLIC "${CMAKE_INSTALL_PREFIX}/include")
target_include_directories(${T} PRIVATE $ENV{EIGEN3_INCLUDE_DIR})
set_property(TARGET ${T} PROPERTY CXX_STANDARD 11)
target_link_libraries(${T} PUBLIC m)
target_link_libraries(${T} PUBLIC stdc++)
install(EXPORT export_${T} FILE overlapConfig.cmake DESTINATION ${CMAKEPREFIX})
add_dependencies(${T} eigen)
