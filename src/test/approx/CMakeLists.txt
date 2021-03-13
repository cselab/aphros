get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} aphros)
add_test_current(NAME "" COMMAND ./test main)

set(T t.${name}.linear)
add_executable(${T} linear.cpp)
target_link_libraries(${T} aphros)
add_test_current(NAME linear COMMAND ./test linear)
