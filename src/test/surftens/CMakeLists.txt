add_test_current(NAME local_8 COMMAND test -c local -b 8 -n 1)
add_test_current(NAME local_16 COMMAND test -c local -b 16 -n 1)
add_test_current(NAME local_32 COMMAND test -c local -b 32 -n 1)

add_test_current(NAME cubismnc_8 COMMAND test -c cubismnc -b 8 -n 2)
add_test_current(NAME cubismnc_16 COMMAND test -c cubismnc -b 16 -n 2)
add_test_current(NAME cubismnc_32 COMMAND test -c cubismnc -b 32 -n 1)

add_test_current(NAME cubismnc_8_vofm COMMAND test -c cubismnc -b 8 -n 2 --vofm)
add_test_current(NAME cubismnc_16_vofm COMMAND test -c cubismnc -b 16 -n 2 --vofm)
add_test_current(NAME cubismnc_32_vofm COMMAND test -c cubismnc -b 32 -n 1 --vofm)

add_test_current(NAME native_8 COMMAND test -c native -b 8 -n 2)
add_test_current(NAME native_16 COMMAND test -c native -b 16 -n 2)
add_test_current(NAME native_32 COMMAND test -c native -b 32 -n 1)
