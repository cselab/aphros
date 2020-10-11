set(T parser)
add_library(${T} parser.cpp)
target_link_libraries(${T} vars filesystem)

set(T vars)
add_library(${T} vars.cpp)

set(T argparse)
add_library(${T} argparse.cpp)
target_link_libraries(${T} vars format)
