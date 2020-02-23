add_test(NAME bc-periodic  COMMAND run periodic)
add_test(NAME bc-slipwall  COMMAND run slipwall)
add_test(NAME bc-surfsymm  COMMAND run surfsymm)
add_test(NAME bc-surfwall  COMMAND run surfwall)
add_test(NAME bc-wall      COMMAND run wall)
add_test(NAME bc-inlet     COMMAND run inlet)
add_test(NAME bc-inletflux COMMAND run inletflux)

add_test(NAME bc-proj-periodic  COMMAND run periodic proj)
add_test(NAME bc-proj-slipwall  COMMAND run slipwall proj)
add_test(NAME bc-proj-surfsymm  COMMAND run surfsymm proj)
add_test(NAME bc-proj-surfwall  COMMAND run surfwall proj)
add_test(NAME bc-proj-wall      COMMAND run wall proj)
add_test(NAME bc-proj-inlet     COMMAND run inlet proj)
add_test(NAME bc-proj-inletflux COMMAND run inletflux proj)

set_tests_properties(
  bc-periodic bc-slipwall bc-surfsymm bc-surfwall
  bc-wall bc-inlet bc-inletflux
  bc-proj-periodic bc-proj-slipwall bc-proj-surfsymm bc-proj-surfwall
  bc-proj-wall bc-proj-inlet bc-proj-inletflux
  PROPERTIES RESOURCE_LOCK bc)
