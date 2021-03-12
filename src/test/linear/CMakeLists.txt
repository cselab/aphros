get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} aphros)

# FIXME must be progated from dependencies
#       but simply linking linear_obj causes a double call of RegisterModule
get_target_property(D linear_obj COMPILE_DEFINITIONS)
target_compile_definitions(${T} PRIVATE ${D})

function(add name)
  add_test_current(NAME ${name} COMMAND ./test ${name})
endfunction()

add(hypre)
add(conjugate)
add(jacobi)
