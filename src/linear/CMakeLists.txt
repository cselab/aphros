set(TARG "hypre")

include_directories(.)
add_library(${TARG} hypre.cpp)

# hypre
set(HYPRE_DIR ${CHPREFIX})
include_directories(${HYPRE_DIR}/include)
target_link_libraries(${TARG} -L${HYPRE_DIR}/lib -lHYPRE -lm)
