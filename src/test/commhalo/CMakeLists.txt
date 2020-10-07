get_filename_component(name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(EXE t.${name})

include_directories(${S})

add_executable(${EXE} main.cpp)

target_link_libraries(${EXE} hydro)
