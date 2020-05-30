add_test(NAME multivof-cub8 COMMAND run 8)
add_test(NAME multivof-cub16 COMMAND run 16)

set_tests_properties(
  multivof-cub8 multivof-cub16
  PROPERTIES RESOURCE_LOCK multivof)
