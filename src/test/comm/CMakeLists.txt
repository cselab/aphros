get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})
include_directories(${S})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} distrsolver approx_eb linear utillinear)

add_test_current(NAME loc COMMAND ap.run ./${EXE} tloc)
add_test_current(NAME cub COMMAND ap.run ./${EXE} tcub)
add_test_current(NAME cubnc COMMAND ap.run ./${EXE} tcubnc)
