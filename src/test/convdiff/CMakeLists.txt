set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distrsolver solver convdiffi)
     
add_test(NAME ${p}-loc COMMAND ap.run ./${EXE} tloc) 
add_test(NAME ${p}-cub COMMAND ap.run ./${EXE} tcub)

set_tests_properties(
  ${p}-loc ${p}-cub
  PROPERTIES RESOURCE_LOCK ${p})
