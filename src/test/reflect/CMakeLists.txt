get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} aphros)
add_test_current(COMMAND ./run)

