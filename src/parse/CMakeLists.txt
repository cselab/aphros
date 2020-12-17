set(T parser)
add_library(${T} parser.cpp)
target_link_libraries(${T} vars filesystem)

set(T vars)
add_library(${T} vars.cpp)

set(T argparse)
add_library(${T} argparse.cpp)
target_link_libraries(${T} vars format)

set(T conf2py)
add_executable(${T} conf2py.cpp)
install(TARGETS ${T} DESTINATION bin)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
target_link_libraries(${T} parser vars format filesystem argparse)

