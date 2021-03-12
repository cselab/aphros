get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(T t.${name})
add_executable(${T} main.cpp)
target_link_libraries(${T} aphros use_dims)

# FIXME must be progated from dependencies
#       but simply linking linear_obj causes a double call of RegisterModule
get_target_property(D linear_obj COMPILE_DEFINITIONS)
target_compile_definitions(${T} PRIVATE ${D})

add_test_current(COMMAND ap.run ./${T})
