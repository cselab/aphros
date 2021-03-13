get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

foreach(DIM RANGE 2 4)
  set(T t.${name}_${DIM}d)
  add_executable_with_objects(${T} main.cpp mesh format logger)
  target_compile_options(${T} PRIVATE -DDIM=${DIM})
  add_test_current(NAME ${DIM}d COMMAND ./test ${DIM}d)
endforeach()
