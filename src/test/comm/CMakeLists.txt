set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distrsolver approx)

add_test(NAME ${p}-loc COMMAND ap.run ./${EXE} tloc)
add_test(NAME ${p}-cub COMMAND ap.run ./${EXE} tcub)
add_test(NAME ${p}-cubnc COMMAND ap.run ./${EXE} tcubnc)

set_tests_properties(
  ${p}-loc ${p}-cub ${p}-cubnc
  PROPERTIES RESOURCE_LOCK ${p})
