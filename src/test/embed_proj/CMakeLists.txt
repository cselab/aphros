set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distrbasic solver proj_eb)

add_test(NAME ${p}-loc COMMAND ap-run ./${EXE} a.conf)
