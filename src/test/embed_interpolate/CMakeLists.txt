get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name})
add_executable(${T} main.cpp)
install(TARGETS ${T} DESTINATION "${APHROS_PREFIX}/bin")
target_link_libraries(${T} aphros)
