cmake_minimum_required(VERSION 3.3.0)
project("aphros")

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/prefix" CACHE PATH "Install path prefix. Defaults to 'prefix' in build directory." FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

include(GNUInstallDirs)
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_INSTALL_BINDIR}")
endif ()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_INSTALL_LIBDIR}")
endif ()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_INSTALL_LIBDIR}")
endif ()

set(SCRIPTS
  applycopyright
  argmax
  argmaxall
  cawk
  copyright
  deriv
  diff
  diffnum
  downsample
  downsampleall
  findsource
  format
  genxmf
  genxmfall
  getcol
  gettraj
  gettrajcol
  gettrajm
  gettrajsh
  getval
  hausdorff
  hist
  histbin
  intcol
  kedr
  line
  makearg
  order
  part
  part0
  plain2vtk
  scale
  spawn
  statcopy
  task
  trajaxes
  vtkcolorall
)
list(TRANSFORM SCRIPTS PREPEND scripts/)

set(SCRIPTS_GEN
  sysinfo
)
list(TRANSFORM SCRIPTS_GEN PREPEND scripts_gen/)

install(PROGRAMS ${SCRIPTS} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
install(PROGRAMS ${SCRIPTS_GEN} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
