cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project("pr")

set(EXE "ex")

aux_source_directory(./ SRC)

include_directories(SRC)

# Add target
add_executable(${EXE} ${SRC})
     

#enable_testing()
include(CTest)

add_test(t1 ex a b)
add_test(t2 ex a a)

