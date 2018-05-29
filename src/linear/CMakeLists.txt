set(TARG "libhypre")

include_directories(.)
add_library(libhypre hypre.cpp)

# hypre
set(HYPRE_DIR ${CHPREFIX})
include_directories(${HYPRE_DIR}/include)
target_link_libraries(${EXE} -lHYPRE -lm)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${HYPRE_DIR}/lib")
