if (USE_OPENCL)
  set(T opencl)
  add_object(${T} opencl.cpp)
  object_link_libraries(${T} mesh use_opencl)

  add_custom_command(
    OUTPUT  kernels.inc
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/kernels.gen"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    DEPENDS kernels.cl
  )
  add_custom_target(kernels_inc DEPENDS kernels.inc)
  set_source_files_properties(main.cpp PROPERTIES OBJECT_DEPENDS kernels.inc)
  add_dependencies(${T} kernels_inc)

  find_package(OpenCL)
  if (OpenCL_FOUND)
    object_link_libraries(${T} OpenCL::OpenCL)
  else()
    find_package(CUDAToolkit QUIET)
    if (CUDAToolkit_FOUND)
      object_link_libraries(${T} CUDA::OpenCL)
    else()
      object_link_libraries(${T} OpenCL)
    endif()
  endif()
endif()
