set(TARG "overlap")

if (OVERLAP)
  set(S "over.cpp")
  find_package(Eigen3 3.0.0 REQUIRED)
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  set(S "overd.cpp")
endif()

add_library(${TARG} ${S})


