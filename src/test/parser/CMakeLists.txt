set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} vars parser)
add_test(NAME ${p} COMMAND ./test_cmp )

set(EXE t.${p}.codeblocks)
add_executable(${EXE} codeblocks.cpp)
add_test(NAME ${p}.codeblocks COMMAND ./run_codeblocks )
