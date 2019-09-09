add_test(NAME surftens-loc8 COMMAND run loc 8)
add_test(NAME surftens-loc16 COMMAND run loc 16)
add_test(NAME surftens-loc32 COMMAND run loc 32 1)
add_test(NAME surftens-cub8 COMMAND run cub 8)
add_test(NAME surftens-cub16 COMMAND run cub 16)
add_test(NAME surftens-cub32 COMMAND run cub 32 1)
add_test(NAME surftens-cubm8 COMMAND run cubm 8)
add_test(NAME surftens-cubm16 COMMAND run cubm 16)
add_test(NAME surftens-cubm32 COMMAND run cubm 32 1)

set_tests_properties(
  surftens-loc8 surftens-loc16 surftens-loc32
  surftens-cub8 surftens-cub16 surftens-cub32
  surftens-cubm8 surftens-cubm16 surftens-cubm32
  PROPERTIES RESOURCE_LOCK surftens)
