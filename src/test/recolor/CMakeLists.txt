set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)

add_test(NAME ${p}-loc COMMAND ch.run ./${EXE} tloc)
add_test(NAME ${p}-cub COMMAND ch.run ./${EXE} tcub)

set_tests_properties(
  ${p}-loc ${p}-cub
  PROPERTIES RESOURCE_LOCK ${p})
