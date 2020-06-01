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

# Scripts

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
install(PROGRAMS ${SCRIPTS} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")


# Select profile
set(GenPath "${CMAKE_CURRENT_SOURCE_DIR}/scripts_gen")
set(ProfilesDir "${GenPath}/profiles")
set(ProfileDefault "default")
file(GLOB ProfileValues RELATIVE "${ProfilesDir}" "${ProfilesDir}/*")
list(PREPEND ProfileValues "${ProfileDefault}")
set(Profile "${ProfileDefault}" CACHE STRING "Profile for tools and wrappers, options are: ${ProfileValues}")
set_property(CACHE Profile PROPERTY STRINGS ${ProfileValues})

# Generated scripts
set(GenOutputDir "${CMAKE_CURRENT_BINARY_DIR}/gen")

if ("${Profile}" STREQUAL "${ProfileDefault}")
  set(GenWorkDir "${GenPath}")
else()
  set(GenWorkDir "${ProfilesDir}/${Profile}")
endif()

set(GenList
  aconf
  base
  conf
  confdeb
  kill
  mpirun
  submit
)
set(GenOutputList "${GenList}")
list(TRANSFORM GenOutputList PREPEND "${GenOutputDir}/")
message("${GenOutputList}")

file(MAKE_DIRECTORY "${GenOutputDir}")

foreach(name ${GenList})
  if (EXISTS "${GenWorkDir}/${name}")
    set(gen "${GenWorkDir}/${name}")
  else()
    set(gen "${GenPath}/${name}")
  endif()
  if (EXISTS "${GenWorkDir}/${name}_inc")
    set(inc "${GenWorkDir}/${name}_inc")
  elseif(EXISTS "${GenPath}/${name}_inc")
    set(inc "${GenPath}/${name}_inc")
  else()
    set(inc)
  endif()
  set(deps "${gen}" "${inc}")
  add_custom_command(OUTPUT "${GenOutputDir}/${name}"
      COMMAND "${gen}" > "${GenOutputDir}/${name}"
      DEPENDS ${deps}
      WORKING_DIRECTORY "${GenWorkDir}"
      COMMENT "Generating script '${name}'")
endforeach()

add_custom_target(scripts_gen ALL DEPENDS "${GenOutputList}")
install(PROGRAMS ${GenOutputList} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
