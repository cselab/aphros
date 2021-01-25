set(T "overlap")

if (NOT DEFINED OVERLAP)
  set(OVERLAP 1)
endif()

if (OVERLAP)
  message("Enable overlap")
  set(S "overlap.cpp")
  include_directories($ENV{EIGEN3_INCLUDE_DIR})
else()
  message("Disable overlap")
  set(S "overlap0.cpp")
endif()

add_library(${T} ${S})
