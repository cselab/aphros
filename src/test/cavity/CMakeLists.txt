add_test(NAME cavity-loc8 COMMAND test_cmp aloc.conf 8)
add_test(NAME cavity-loc16 COMMAND test_cmp aloc.conf 16)
add_test(NAME cavity-loc32 COMMAND test_cmp aloc.conf 32 1)
add_test(NAME cavity-cub8 COMMAND test_cmp acub.conf 8)
add_test(NAME cavity-cub16 COMMAND test_cmp acub.conf 16)
add_test(NAME cavity-cub32 COMMAND test_cmp acub.conf 32 1)

set_tests_properties(
  cavity-loc8 cavity-loc16 cavity-loc32
  cavity-cub8 cavity-cub16 cavity-cub32
  PROPERTIES RESOURCE_LOCK cavity)
