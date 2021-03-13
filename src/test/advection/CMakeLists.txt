get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
install(TARGETS ${T} DESTINATION "${APHROS_PREFIX}/bin")
target_link_libraries(${T} aphros)


foreach(comm local cubismnc native)
  foreach(block 8 16 32)
    add_test_current(
        NAME 2d-${comm}${block}
        COMMAND ./test 2d -b ${block} -c ${comm})
  endforeach()
  add_test_current(
      NAME 3d-${comm}${block}
      COMMAND ./test 3d -c ${comm})
endforeach()
