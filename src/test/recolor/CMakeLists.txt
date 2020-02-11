set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distrbasic overlap utilvof)

add_test(NAME ${p} COMMAND ap.run ./${EXE})
