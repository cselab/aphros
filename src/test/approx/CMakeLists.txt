set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} suspender approx timer approx_eb)
add_test(NAME ${p} COMMAND ./run)

set(EXE t.${p}.linear)
add_executable(${EXE} linear.cpp)
target_link_libraries(${EXE} suspender approx_eb)
add_test(NAME ${p}.linear COMMAND ./run_linear)
