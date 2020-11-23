get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(EXE t.${name})
include_directories(${S})
add_executable(${EXE} main.cpp)
target_link_libraries(${EXE} aphros)

add_test_current(NAME local COMMAND ap.mpirun -n 1  ./${EXE} --extra "set string backend local")
add_test_current(NAME cubismnc COMMAND ap.mpirun -n 2  ./${EXE} --extra "set string backend cubismnc")
