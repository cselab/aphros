cmake_minimum_required(VERSION 3.3.0)
project("aphros")

set(NAMEPREFIX ap.)

if (DEFINED ENV{APHROS_PREFIX})
  set(CMAKE_INSTALL_PREFIX $ENV{APHROS_PREFIX}
      CACHE PATH "Install path prefix. Set from APHROS_PREFIX" FORCE)
elseif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/prefix"
      CACHE PATH "Install path prefix. Defaults to 'prefix' in build directory." FORCE)
endif ()

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
set(
    Profile "${ProfileDefault}"
    CACHE STRING "Profile for tools and wrappers, options are: ${ProfileValues}")
set_property(CACHE Profile PROPERTY STRINGS ${ProfileValues})

# Scripts
set(ScriptsDir "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
set(ScriptsProfileDir "${ScriptsDir}/profiles/${Profile}")
set(ScriptsList
    aconf
    applycopyright
    argmax
    argmaxall
    base
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
    kedr
    kill
    line
    makearg
    makesim
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
set(ScriptsListNoExec
    job
    job0
    sim_base.makefile
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

foreach(name ${ScriptsListNoExec})
  if (EXISTS "${ScriptsProfileDir}/${name}")
    set(path "${ScriptsProfileDir}/${name}")
  else()
    set(path "${ScriptsDir}/${name}")
  endif()

  install(
      FILES "${path}"
      RENAME ${NAMEPREFIX}${name}
      DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
endforeach()


install(
    DIRECTORY
      "${ScriptsDir}/aphros"
      "${ScriptsDir}/paratools"
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
    create_a_conf
    create_base_conf
    setenv
    submit
  )
set(GenListNoExec
    setenv
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
  add_custom_command(
      OUTPUT "${output}"
      COMMAND ${CMAKE_COMMAND} -E env PREFIX=${CMAKE_INSTALL_PREFIX} "${gen}" > "${output}"
      DEPENDS ${deps}
      WORKING_DIRECTORY "${GenWorkDir}"
      COMMENT "Generating script '${name}'")
  if ("${name}" IN_LIST GenListNoExec)
    install(
        FILES "${output}"
        RENAME "${NAMEPREFIX}${name}"
        DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  else()
    install(
        PROGRAMS "${output}"
        RENAME "${NAMEPREFIX}${name}"
        DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
  endif()
endforeach()
add_custom_target(scripts_gen ALL DEPENDS "${GenOutputList}")

include(ExternalProject)

# Hypre
set(hypre_tar "${CMAKE_CURRENT_SOURCE_DIR}/hypre/hypre.tar.gz")
ExternalProject_Add(hypre
    URL ${hypre_tar}
    CMAKE_ARGS
      -DHYPRE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DHYPRE_SHARED=ON
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/hypre
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    SOURCE_SUBDIR src
  )
