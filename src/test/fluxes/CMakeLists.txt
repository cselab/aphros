function(add name)
  add_test_current(NAME ${name} COMMAND ./test ${name})
endfunction()

add(simple_imp_gravity)
add(simple_exp_gravity)
add(proj_imp_gravity)
add(proj_exp_gravity)
add(proj_imp_pressure)
