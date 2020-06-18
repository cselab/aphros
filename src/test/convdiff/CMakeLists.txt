get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} distrsolver solver convdiffi)

add_test_current(NAME loc COMMAND ap.run ./${EXE} tloc)
add_test_current(NAME cub COMMAND ap.run ./${EXE} tcub)
