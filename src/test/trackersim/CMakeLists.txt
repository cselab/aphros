get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test(NAME ${p}-zero COMMAND ./test zero)
add_test(NAME ${p}-uniform COMMAND ./test uniform)
add_test(NAME ${p}-gravity COMMAND ./test gravity)
add_test(NAME ${p}-grad COMMAND ./test grad)

set_tests_properties(
  ${p}-zero ${p}-uniform ${p}-gravity ${p}-grad
  PROPERTIES RESOURCE_LOCK tracker)
