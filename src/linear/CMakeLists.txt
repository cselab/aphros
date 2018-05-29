set(TARG "hypre")

#include_directories(.)
add_library(${TARG} hypre.cpp)

# hypre
set(HYPRE_DIR ${CHPREFIX})
target_include_directories(${TARG} PUBLIC ${HYPRE_DIR}/include)
target_link_libraries(${TARG} -L${HYPRE_DIR}/lib -lHYPRE -lm)
