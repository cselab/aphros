get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(EXE t.${name}.getmap)
add_executable(${EXE} getmap.cpp)
target_link_libraries(${EXE} aphros)
add_test_current(NAME getmap COMMAND ./test getmap)

set(EXE t.${name}.parse)
add_executable(${EXE} parse.cpp)
target_link_libraries(${EXE} aphros)
add_test_current(NAME parse COMMAND ./test parse)
