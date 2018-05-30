set(CSD "${CMAKE_CURRENT_SOURCE_DIR}")
get_filename_component(p ${CSD} NAME)

set(EXE t.${p})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} distr)
     
add_test(NAME ${p}-loc COMMAND mpirun -n 2 $<TARGET_FILE:${EXE}> loc.conf WORKING_DIRECTORY "${CSD}")
add_test(NAME ${p}-cub COMMAND mpirun -n 2 $<TARGET_FILE:${EXE}> cub.conf WORKING_DIRECTORY "${CSD}")
