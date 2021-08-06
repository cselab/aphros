get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})

add_executable(${T} main.cpp)
target_link_libraries(${T} aphros overlap)
add_test_current(NAME "" COMMAND ${T})

add_executable(plane plane.cpp)

set(T levelset)
add_executable(${T} levelset.cpp)
target_link_libraries(${T} use_mpi)
add_test_current(NAME ${T} COMMAND run_levelset)
