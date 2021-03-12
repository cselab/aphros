get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp ../../inside/main.c)
target_link_libraries(${T} sysinfo approx_eb parser primlist inside)
target_compile_options(${T} PRIVATE -O3 -march=native -ffast-math)

add_test_current(COMMAND ${T})

