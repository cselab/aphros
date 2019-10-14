set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distrbasic)

add_test(NAME ${p} COMMAND ch.run ./${EXE})

set_tests_properties(
  ${p}
  PROPERTIES RESOURCE_LOCK ${p})
