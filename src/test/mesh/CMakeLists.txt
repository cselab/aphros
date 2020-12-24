get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name}_2d)
add_executable(${T} main.cpp)
target_compile_options(${T} PRIVATE -DDIM=2)
target_link_libraries(${T} aphros)

set(T t.${name}_3d)
add_executable(${T} main.cpp)
target_compile_options(${T} PRIVATE -DDIM=3)
target_link_libraries(${T} aphros)

set(T t.${name}_4d)
add_executable(${T} main.cpp)
target_compile_options(${T} PRIVATE -DDIM=4)
target_link_libraries(${T} aphros)

add_test_current(NAME 2d COMMAND ./test 2d)
add_test_current(NAME 3d COMMAND ./test 3d)
add_test_current(NAME 4d COMMAND ./test 4d)
