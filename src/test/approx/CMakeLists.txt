get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(EXE t.${name})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} hydro)
add_test_current(NAME "" COMMAND ./test main)

set(EXE t.${name}.linear)
add_executable(${EXE} linear.cpp)
target_link_libraries(${EXE} hydro)
add_test_current(NAME linear COMMAND ./test linear)
