add_test(NAME contang_2d COMMAND run 2d)
add_test(NAME contang_3d COMMAND run 3d)
add_test(NAME contang_embed_2d COMMAND run embed_2d)
add_test(NAME contang_embed_3d COMMAND run embed_3d)

set_tests_properties(
  contang_2d contang_3d contang_embed_2d contang_embed_3d
  PROPERTIES RESOURCE_LOCK contang)
