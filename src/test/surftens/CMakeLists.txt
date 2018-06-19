add_test(NAME surftens-loc8 COMMAND test_cmp aloc.conf 8)
add_test(NAME surftens-loc16 COMMAND test_cmp aloc.conf 16)
add_test(NAME surftens-loc32 COMMAND test_cmp aloc.conf 32 1)
add_test(NAME surftens-cub8 COMMAND test_cmp acub.conf 8)
add_test(NAME surftens-cub16 COMMAND test_cmp acub.conf 16)
add_test(NAME surftens-cub32 COMMAND test_cmp acub.conf 32 1)

