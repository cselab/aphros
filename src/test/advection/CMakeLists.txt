get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})
add_executable(${EXE} main.cpp)
install(TARGETS ${EXE} DESTINATION "${APHROS_PREFIX}/bin")
target_link_libraries(${EXE} hydro)


foreach(comm local cubismnc cubism)
  foreach(block 8 16 32)
    add_test_current(
        NAME 2d-${comm}${block}
        COMMAND ./test 2d -b ${block} -c ${comm})
  endforeach()
  add_test_current(
      NAME 3d-${comm}${block}
      COMMAND ./test 3d -c ${comm})
endforeach()
