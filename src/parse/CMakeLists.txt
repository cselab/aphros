add_library(parser parser.cpp)
target_link_libraries(parser vars filesystem)

add_library(vars vars.cpp)
