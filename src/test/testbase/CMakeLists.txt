get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test(NAME ${name}.b0 COMMAND ./test case0)
add_test(NAME ${name}.b1 COMMAND ./test case1)
