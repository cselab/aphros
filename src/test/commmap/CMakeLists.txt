get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (USE_MPI)
  set(T t.${name})
  add_executable(${T} main.cpp)
  target_link_libraries(${T} aphros openmp)
  add_test_current(COMMAND ./test main)
  # With more threads than ranks, the number of subdomains
  # (ranks / threads) is zero, and the test aborts.
  set_tests_properties(${name} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)
endif()

set(T t.${name}.manager)
add_executable(${T} manager.cpp)
target_link_libraries(${T} aphros)
if (USE_MPI)
  add_test_current(NAME manager COMMAND ./test manager)
else()
  add_test_current(NAME manager COMMAND ./test manager_seq)
endif()

set(T t.${name}.rank)
add_executable(${T} rank.cpp)
target_link_libraries(${T} aphros)
if (USE_BACKEND_CUBISM)
  add_test_current(NAME rank-cubism COMMAND ./test rank_cubism)
endif()
add_test_current(NAME rank-native COMMAND ./test rank_native)
