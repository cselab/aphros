set(T parser)
add_object(${T} parser.cpp codeblocks.cpp)
object_link_libraries(${T} vars filesystem)

set(T vars)
add_object(${T} vars.cpp)
object_link_libraries(${T} logger)

set(T evalexpr)
add_object(${T} evalexpr.cpp)
object_link_libraries(${T} logger format)

set(T parse_template)
add_object(${T} template.cpp)

set(T argparse)
add_object(${T} argparse.cpp)
object_link_libraries(${T} vars format)
