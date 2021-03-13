get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name}.getmap)
add_executable(${T} getmap.cpp)
target_link_libraries(${T} aphros)
add_test_current(NAME getmap COMMAND ./test getmap)

set(T t.${name}.parse)
add_executable(${T} parse.cpp)
target_link_libraries(${T} aphros)
add_test_current(NAME parse COMMAND ./test parse)
