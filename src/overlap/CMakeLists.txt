set(TARG "overlap")

if (OVERLAP)
  set(S "over.cpp")
  include_directories($ENV{EIGEN3_INCLUDE_DIR})
else()
  set(S "overd.cpp")
endif()

add_library(${TARG} ${S})
