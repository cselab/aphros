get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name})
add_executable_with_objects(${T} main.cpp argparse format parser)
add_test_current(NAME main COMMAND ./test main)

set(T t.${name}.codeblocks)
add_executable_with_objects(${T} codeblocks.cpp argparse format parser)
add_test_current(NAME codeblocks COMMAND ./test codeblocks)
