get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${p})
add_executable_with_objects(${T} main.cpp mesh format logger)
add_test_current(COMMAND test)

