get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} hydro)
add_test_current(NAME main COMMAND ./test main)

set(T t.${name}.codeblocks)
add_executable(${T} codeblocks.cpp)
add_test_current(NAME codeblocks COMMAND ./test codeblocks)
