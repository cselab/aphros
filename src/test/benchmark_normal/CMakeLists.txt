get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} suspender format timer sysinfo mesh allocator)
target_compile_options(${T} PRIVATE -march=native)
add_test_current(COMMAND ${T})

