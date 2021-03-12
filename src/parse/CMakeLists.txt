set(T parser)
aphros_add_object(T ${T} parser.cpp codeblocks.cpp)
target_link_libraries(${T} vars filesystem)

set(T vars)
aphros_add_object(T ${T} vars.cpp)
target_link_libraries(${T} logger)

set(T parse_template)
aphros_add_object(T ${T} template.cpp)

set(T argparse)
aphros_add_object(T ${T} argparse.cpp)
target_link_libraries(${T} vars format)

set(T conf2py)
add_executable(${T} conf2py.cpp)
install(TARGETS ${T} DESTINATION bin)
set_target_properties(${T} PROPERTIES PREFIX "ap.")
target_link_libraries(${T} parser vars format filesystem argparse logger)
