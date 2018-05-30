set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)
     
add_test(NAME ${p}-loc COMMAND mpirun -n 2 ${EXE} loc.conf)
add_test(NAME ${p}-cub COMMAND mpirun -n 2 ${EXE} cub.conf)
