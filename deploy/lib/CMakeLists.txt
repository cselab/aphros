add_subdirectory(table)
add_subdirectory(vof)
add_subdirectory(vofi)
add_subdirectory(vtk)
add_subdirectory(vtkcolor)

if (USE_EIGEN)
  add_subdirectory(overlap)
endif()

if (USE_GSL)
  add_subdirectory(chmath)
  add_subdirectory(csv)
  add_subdirectory(split)
endif()

if (USE_HDF5)
  add_subdirectory(h5)
endif()

