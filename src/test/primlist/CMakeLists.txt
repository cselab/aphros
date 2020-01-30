set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE0 t.${p}.getmap)
set(EXE1 t.${p}.parse)

add_executable(${EXE0} getmap.cpp)
add_executable(${EXE1} parse.cpp)

target_link_libraries(${EXE0} vars parser)
target_link_libraries(${EXE1} vars parser)

add_test(NAME ${p}_getmap COMMAND ./run_getmap)
add_test(NAME ${p}_parse COMMAND ./run_parse)
