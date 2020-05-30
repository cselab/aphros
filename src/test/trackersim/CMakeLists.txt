add_test(NAME tracker-zero COMMAND run zero)
add_test(NAME tracker-uniform COMMAND run uniform)
add_test(NAME tracker-gravity COMMAND run gravity)
add_test(NAME tracker-grad COMMAND run grad)

set_tests_properties(
  tracker-zero tracker-uniform tracker-gravity tracker-grad
  PROPERTIES RESOURCE_LOCK tracker)
