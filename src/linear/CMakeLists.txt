set(T "linear")
add_library(${T} linear.cpp)
target_link_libraries(${T})
target_compile_definitions(${T} PUBLIC _USE_HYPRE_=$<BOOL:${USE_HYPRE}>)

if (USE_HYPRE)
  set(T "hypre")
  add_library(${T} hypre.cpp)
  target_link_libraries(${T} hypreext)

  set(T "hypresub")
  add_library(${T} hypresub.cpp)
  target_link_libraries(${T} hypre subcomm histogram)

  set(T "linear_hypre")
  add_library(${T} linear_hypre.cpp)
  target_link_libraries(${T} hypre)
endif()
