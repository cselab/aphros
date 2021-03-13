get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable_with_objects(${T} main.cpp sysinfo approx_eb parser primlist)
target_link_libraries(${T} inside openmp)
target_compile_options(${T} PRIVATE -O3 -march=native -ffast-math)

add_test_current(COMMAND ${T})

