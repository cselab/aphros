set(T "linear")
add_object(${T} linear.cpp)
object_link_libraries(${T} use_mpi use_dims)
object_compile_definitions(${T} PUBLIC _USE_HYPRE_=$<BOOL:${USE_HYPRE}>)
object_compile_definitions(${T} PUBLIC _USE_AMGX_=$<BOOL:${USE_AMGX}>)

if (USE_HYPRE)
  set(T "hypre")
  add_object(${T} hypre.cpp)
  object_link_libraries(${T} hypreext use_mpi)

  set(T "hypresub")
  add_object(${T} hypresub.cpp)
  object_link_libraries(${T} hypre subcomm histogram)

  set(T "linear_hypre")
  add_object(${T} linear_hypre.cpp)
  object_link_libraries(${T} hypre)
endif()

if (USE_AMGX)
  set(T "linear_amgx")
  add_object(${T} linear_amgx.cpp)
  target_include_directories(${T} PRIVATE ${APHROS_PREFIX}/include)
  object_link_libraries(${T} amgxsh use_mpi)
  if (CUDAToolkit_FOUND)
    object_link_libraries(${T} CUDA::cudart_static)
  endif()
endif()

if (USE_OPENCL)
  set(T "linear_conjugate_cl")
  add_object(${T} conjugate_cl.cpp)
  object_link_libraries(${T} use_mpi use_dims use_opencl)
  add_custom_command(
    OUTPUT  conjugate_cl.inc
    COMMAND ./kernel.gen conjugate_cl.cl conjugate_cl.inc
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS conjugate_cl.cl
  )
  add_custom_target(conjugate_cl_inc DEPENDS conjugate_cl.inc)
  set_source_files_properties(cojnugate_cl.cpp
      PROPERTIES OBJECT_DEPENDS conjugate_cl.inc)
  add_dependencies(${T} conjugate_cl_inc)

  find_package(OpenCL)
  if (OpenCL_FOUND)
    object_link_libraries(${T} PRIVATE OpenCL::OpenCL)
  else()
    find_package(CUDAToolkit QUIET)
    if (CUDAToolkit_FOUND)
      object_link_libraries(${T} PRIVATE CUDA::OpenCL)
    else()
      object_link_libraries(${T} PRIVATE OpenCL)
    endif()
  endif()
endif()
