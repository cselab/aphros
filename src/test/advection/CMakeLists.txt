set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")

target_link_libraries(${EXE} distr solver vof tvd)

add_test(NAME adv-loc8 COMMAND runshort tloc 8)
add_test(NAME adv-loc16 COMMAND runshort tloc 16)
add_test(NAME adv-loc32 COMMAND runshort tloc 32 1)
add_test(NAME adv-cub8 COMMAND runshort tcub 8)
add_test(NAME adv-cub16 COMMAND runshort tcub 16)
add_test(NAME adv-cub32 COMMAND runshort tcub 32 1)

