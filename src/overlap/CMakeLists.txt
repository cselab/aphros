if (USE_OVERLAP)
  set(T "overlap")
  add_library(${T} overlap.cpp)
  target_include_directories(${T} PRIVATE $ENV{EIGEN3_INCLUDE_DIR})
endif()

