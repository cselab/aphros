if (USE_OVERLAP)
  set(T "overlap")
  add_object(${T} overlap.cpp)
  target_include_directories(${T}_obj PRIVATE $ENV{EIGEN3_INCLUDE_DIR})
endif()

