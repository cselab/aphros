cmake_minimum_required(VERSION 3.3.0)

project("mfer")

set(EXE "mfer")

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wno-deprecated-copy")

add_executable(${EXE} main.cpp)
set_target_properties(${EXE} PROPERTIES PREFIX "ch.")

# git revision
set(P "${PROJECT_SOURCE_DIR}/..")
set(C "${P}/src/util/gitgen.cpp")
# create on configure
if (NOT EXISTS ${C})
  execute_process(
    COMMAND "${P}/tool/gitrev" "${C}"
    WORKING_DIRECTORY "${P}"
  )
endif()
# create on build
add_custom_target(gitrev
  COMMAND "${P}/tool/gitrev" "${C}"
  WORKING_DIRECTORY "${P}"
)
add_dependencies(${EXE} gitrev)

# install prefix
set(CHPREFIX $ENV{CHPREFIX})
install(TARGETS ${EXE} DESTINATION "${CHPREFIX}/bin")

include_directories(.)

add_subdirectory(linear)
add_subdirectory(distr)
add_subdirectory(util)
add_subdirectory(solver)
add_subdirectory(parse)
add_subdirectory(dump)
add_subdirectory(overlap)
add_subdirectory(young)

target_link_libraries(${EXE} 
  distr solver vof tvd convdiffvi simple overlap young utilfluid events)

add_subdirectory(test)
enable_testing()

if(LIB)
	set(HYPRE_DIR ${CHPREFIX})
	set(T hydro)
	include(header.cmake)
	include(lib.cmake)
	target_compile_definitions(${T} PUBLIC _ALIGNBYTES_=16)
	target_include_directories(${T} PUBLIC ${HYPRE_DIR}/include)
	target_link_libraries(${T} -L${HYPRE_DIR}/lib -lHYPRE -lm)
	target_include_directories(${T} PUBLIC ${HDF5_INCLUDE_DIRS})
	target_link_libraries(${T} ${HDF5_LIBRARIES})
	target_compile_definitions(${T} PUBLIC _USE_HDF_)
	install(TARGETS ${T} ARCHIVE DESTINATION ${CHPREFIX}/lib)
endif(LIB)
