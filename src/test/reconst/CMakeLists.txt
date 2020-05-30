set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} overlap)

add_test(NAME ${p} COMMAND ${EXE})

add_executable(plane plane.cpp)

add_executable(levelset levelset.cpp)
add_test(NAME ${p}.levelset COMMAND run_levelset)
