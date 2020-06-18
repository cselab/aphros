get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} suspender sysinfo timer)
add_test_current(COMMAND ${EXE})

