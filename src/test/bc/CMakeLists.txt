add_test(NAME bc-periodic  COMMAND run periodic)
add_test(NAME bc-slipwall  COMMAND run slipwall)
add_test(NAME bc-surfsymm  COMMAND run surfsymm)
add_test(NAME bc-surfwall  COMMAND run surfwall)
add_test(NAME bc-wall      COMMAND run wall)
add_test(NAME bc-inlet     COMMAND run inlet)
add_test(NAME bc-inletflux COMMAND run inletflux)

set_tests_properties(
  bc-periodic bc-slipwall bc-surfsymm bc-surfwall
  bc-wall bc-inlet bc-inletflux
  PROPERTIES RESOURCE_LOCK bc)
