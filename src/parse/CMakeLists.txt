set(T parser)
add_library(${T} parser.cpp codeblocks.cpp)
target_link_libraries(${T} vars filesystem utilsystem)

set(T vars)
add_library(${T} vars.cpp)

set(T parse_template)
add_library(${T} template.cpp)

set(T argparse)
add_library(${T} argparse.cpp)
target_link_libraries(${T} vars format)

set(T conf2py)
add_executable(${T} conf2py.cpp)
install(TARGETS ${T} DESTINATION bin)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
target_link_libraries(${T} parser vars format filesystem argparse logger utilsystem)

