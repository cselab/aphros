get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable_with_objects(${T} main.cpp format timer sysinfo mesh allocator openmp)
target_link_libraries(${T} use_avx)
add_test_current(COMMAND ${T})

