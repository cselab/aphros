add_library(parser parser.cpp)
target_link_libraries(parser vars)

add_library(vars vars.cpp)
