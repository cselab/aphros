cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(EXE "main")

set(H "../../hydro")

include_directories(${H})

aux_source_directory(${H} SRC)

add_executable(${EXE} main.cpp ${SRC})
     
set_property(TARGET ${EXE} PROPERTY CXX_STANDARD 11)

enable_testing()

add_test(main ${EXE} a b)

