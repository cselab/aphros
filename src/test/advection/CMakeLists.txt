set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)

add_test(NAME adv-loc8 COMMAND test_cmp al8.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME adv-loc16 COMMAND test_cmp al16.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME adv-loc32 COMMAND test_cmp al32.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME adv-cub8 COMMAND test_cmp alc8.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME adv-cub16 COMMAND test_cmp alc16.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME adv-cub32 COMMAND test_cmp alc32.conf 1 WORKING_DIRECTORY "${CSD}")

