get_filename_component(p ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)
     
enable_testing()

add_test(${p}-loc mpirun -n 2 ${EXE} loc.conf)
add_test(${p}-cub mpirun -n 2 ${EXE} cub.conf)
