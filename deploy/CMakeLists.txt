cmake_minimum_required(VERSION 3.3.0)
project("aphros")

set(NAMEPREFIX ap.)

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

# Select profile
# FIXME: combine scripts/profiles and scripts_gen/profiles
#        now only looking in "scripts_gen/profiles"
set(GenDir "${CMAKE_CURRENT_SOURCE_DIR}/scripts_gen")
set(ProfilesDir "${GenDir}/profiles")
set(ProfileDefault "default")
file(GLOB ProfileValues RELATIVE "${ProfilesDir}" "${ProfilesDir}/*")
list(APPEND ProfileValues "${ProfileDefault}")
set(Profile "${ProfileDefault}" CACHE STRING "Profile for tools and wrappers, options are: ${ProfileValues}")
set_property(CACHE Profile PROPERTY STRINGS ${ProfileValues})

# Scripts
set(ScriptsDir "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
set(ScriptsProfileDir "${ScriptsDir}/profiles/${Profile}")
set(ScriptsList
  applycopyright
  argmax
  argmaxall
  cawk
  conf
  confdeb
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
  job
  job0
  kedr
  kill
  line
  makearg
  mpirun
  order
  part
  part0
  plain2vtk
  run
  run0
  scale
  spawn
  statcopy
  submit0
  sysinfo
  task
  trajaxes
  vtkcolorall
  writearg
)

foreach(name ${ScriptsList})
  if (EXISTS "${ScriptsProfileDir}/${name}")
    set(path "${ScriptsProfileDir}/${name}")
  else()
    set(path "${ScriptsDir}/${name}")
  endif()
  install(
    PROGRAMS "${path}"
    RENAME ${NAMEPREFIX}${name}
    DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endforeach()

install(
  DIRECTORY "${ScriptsDir}/aphros" "${ScriptsDir}/paratools"
  DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
  USE_SOURCE_PERMISSIONS
  PATTERN __pycache__ EXCLUDE)

# Generated scripts
set(GenProfileDir "${GenDir}/profiles/${Profile}")
set(GenOutputDir "${CMAKE_CURRENT_BINARY_DIR}/scripts_gen")
if ("${Profile}" STREQUAL "${ProfileDefault}")
  set(GenWorkDir "${GenDir}")
else()
  set(GenWorkDir "${ProfilesDir}/${Profile}")
endif()

set(GenList
  aconf
  base
  submit
)
file(MAKE_DIRECTORY "${GenOutputDir}")
set(GenOutputList)
foreach(name ${GenList})
  if (EXISTS "${GenWorkDir}/${name}")
    set(gen "${GenWorkDir}/${name}")
  else()
    set(gen "${GenDir}/${name}")
  endif()
  if (EXISTS "${GenWorkDir}/${name}_inc")
    set(inc "${GenWorkDir}/${name}_inc")
  elseif(EXISTS "${GenDir}/${name}_inc")
    set(inc "${GenDir}/${name}_inc")
  else()
    set(inc)
  endif()
  set(deps "${gen}" "${inc}")
  set(output "${GenOutputDir}/${name}")
  list(APPEND GenOutputList "${output}")
  add_custom_command(OUTPUT "${output}"
      COMMAND "${gen}" > "${output}"
      DEPENDS ${deps}
      WORKING_DIRECTORY "${GenWorkDir}"
      COMMENT "Generating script '${name}'")
  install(
    PROGRAMS "${output}"
    RENAME "${NAMEPREFIX}${name}"
    DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endforeach()
add_custom_target(scripts_gen ALL DEPENDS "${GenOutputList}")
