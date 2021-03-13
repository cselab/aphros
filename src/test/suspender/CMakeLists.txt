get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable_with_objects(${T} main.cpp suspender logger)
add_test_current(COMMAND ${T})

