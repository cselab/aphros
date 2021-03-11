get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} argparse format parser filesystem logger)
add_test_current(NAME main COMMAND ./test main)

set(T t.${name}.codeblocks)
add_executable(${T} codeblocks.cpp)
target_link_libraries(${T} logger format parser filesystem)
add_test_current(NAME codeblocks COMMAND ./test codeblocks)
