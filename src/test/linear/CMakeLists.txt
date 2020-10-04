get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} distrbasic subcomm approx_eb linear)

function(add name)
  add_test_current(NAME ${name} COMMAND ./test ${name})
endfunction()

add(hypre)
add(conjugate)
add(jacobi)
