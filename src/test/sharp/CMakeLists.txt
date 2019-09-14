add_test(NAME sharp COMMAND run cub)
add_test(NAME sharpm COMMAND run cubm)

set_tests_properties(
  sharp sharpm
  PROPERTIES RESOURCE_LOCK sharp)
