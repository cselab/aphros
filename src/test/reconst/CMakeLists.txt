get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})

add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} overlap)
add_test_current(NAME "" COMMAND ${EXE})

add_executable(plane plane.cpp)

add_executable(levelset levelset.cpp)
add_test_current(NAME levelset COMMAND run_levelset)
