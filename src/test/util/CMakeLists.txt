get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp $<TARGET_OBJECTS:utildistr> $<TARGET_OBJECTS:format>)
add_test_current(COMMAND ./test)
