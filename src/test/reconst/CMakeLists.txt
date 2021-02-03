get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})

add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} aphros)
add_test_current(NAME "" COMMAND ${EXE})

add_executable(plane plane.cpp)

set(T levelset)
add_executable(${T} levelset.cpp)
target_link_libraries(${T} use_mpi)
add_test_current(NAME ${T} COMMAND run_levelset)
