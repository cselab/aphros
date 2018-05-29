get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)

enable_testing()

add_test(NAME adv-loc8 WORKING_DIRECTORY . COMMAND test_cmp al8.conf)
add_test(adv-loc16 test_cmp al16.conf)
add_test(adv-loc32 test_cmp al32.conf)
add_test(adv-cub8 test_cmp alc8.conf)
add_test(adv-cub16 test_cmp alc16.conf)
add_test(adv-cub32 test_cmp alc32.conf 1)

