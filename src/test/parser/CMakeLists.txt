get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${p})
add_executable(${T} main.cpp)
target_link_libraries(${T} vars parser)
add_test(NAME ${p}.main COMMAND ./test main)

set(T t.${p}.codeblocks)
add_executable(${T} codeblocks.cpp)
add_test(NAME ${p}.codeblocks COMMAND ./test codeblocks)
