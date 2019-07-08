add_test(NAME surftens-loc8 COMMAND test_cmp aloc.conf 8)
add_test(NAME surftens-loc16 COMMAND test_cmp aloc.conf 16)
add_test(NAME surftens-loc32 COMMAND test_cmp aloc.conf 32 1)
add_test(NAME surftens-cub8 COMMAND test_cmp acub.conf 8)
add_test(NAME surftens-cub16 COMMAND test_cmp acub.conf 16)
add_test(NAME surftens-cub32 COMMAND test_cmp acub.conf 32 1)

set_tests_properties(
  surftens-loc8 surftens-loc16 surftens-loc32
  surftens-cub8 surftens-cub16 surftens-cub32
  PROPERTIES RESOURCE_LOCK surftens)
