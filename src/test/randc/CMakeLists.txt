add_test(NAME randc-dim2  COMMAND run dim=2)
add_test(NAME randc-dim3  COMMAND run dim=3)

set_tests_properties(
  randc-dim2 randc-dim3
  PROPERTIES RESOURCE_LOCK randc)
