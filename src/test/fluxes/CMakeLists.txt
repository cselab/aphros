function(add a0 a1 a2)
  set(name ${a0}-${a1}-${a2})
  add_test_current(NAME ${name} COMMAND run ${a0} ${a1} ${a2} ${ARGN})
endfunction()

add(simple imp gravity)
add(simple exp gravity)
add(proj imp gravity)
add(proj exp gravity)
