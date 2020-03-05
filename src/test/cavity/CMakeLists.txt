add_test(NAME cavity-loc8 COMMAND run aloc.conf 8)
add_test(NAME cavity-loc16 COMMAND run aloc.conf 16)
add_test(NAME cavity-loc32 COMMAND run aloc.conf 32 1)
add_test(NAME cavity-cub8 COMMAND run acub.conf 8)
add_test(NAME cavity-cub16 COMMAND run acub.conf 16)
add_test(NAME cavity-cub32 COMMAND run acub.conf 32 1)

set_tests_properties(
  cavity-loc8 cavity-loc16 cavity-loc32
  cavity-cub8 cavity-cub16 cavity-cub32
  PROPERTIES RESOURCE_LOCK cavity)
