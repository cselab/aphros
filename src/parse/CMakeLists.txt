add_library(interp interp.cpp)
target_link_libraries(interp vars)

add_library(vars vars.cpp)
