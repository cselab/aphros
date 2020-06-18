get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test(NAME ${p} COMMAND ./test)
