get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)

target_link_libraries(${T} mesh suspender timer sysinfo approx_eb vars parser
  format filesystem approx_eb approx primlist)
target_compile_options(${T} PRIVATE -O3 -march=native)

add_test_current(COMMAND ${T})

