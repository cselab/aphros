get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(T t.${name}_util)
add_executable(${T} dump_util.cpp)
target_link_libraries(${T} aphros)
add_test_current(NAME dump_util COMMAND ./test util)

set(T t.${name}_diff)
add_executable(${T} dump_diff.cpp)
target_link_libraries(${T} aphros)

set(T t.${name}_gen)
add_executable(${T} dump_gen.cpp)
target_link_libraries(${T} aphros)

set(T t.${name}_meta)
add_executable(${T} dump_meta.cpp)
target_link_libraries(${T} aphros)

if (USE_MPI)
  add_test_current(NAME dump_raw COMMAND ./test raw)
  add_test_current(NAME dump_h5 COMMAND ./test h5)
else()
  add_test_current(NAME dump_raw COMMAND ./test raw_seq)
endif()
add_test_current(NAME dump_meta COMMAND ./test meta)
