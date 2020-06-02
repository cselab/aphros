set(T march)
add_library(${T} SHARED
  main.c
  march.h
  )
set_property(TARGET ${T} PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(${T} PROPERTIES PUBLIC_HEADER "march.h")
install(TARGETS ${T} 
    EXPORT export DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
target_include_directories(${T} PUBLIC "${CMAKE_INSTALL_PREFIX}/include")
target_link_libraries(${T} PUBLIC m)
install(EXPORT export FILE marchConfig.cmake DESTINATION ${CMAKEPREFIX})
