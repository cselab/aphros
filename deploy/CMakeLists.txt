cmake_minimum_required(VERSION 3.10.2)
project("aphros-deploy")

option(USE_FPZIP "Install fpzip. Used for compression of communication buffers (experimental)" OFF)
option(USE_HYPRE "Install Hypre library providing linear solvers" ON)
option(USE_AMGX "Install AMGX, linear solver on GPUs. Requires CUDA Toolkit" OFF)
option(USE_HDF5 "Install utilities that require HDF5: h5" OFF)

set(NAMEPREFIX ap.)

set(CMAKE_MACOSX_RPATH OFF)
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY ON)

# Disallow in-source build
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR
    "Aphros requires an out-of-source build. Run CMake from a separate build directory.")
endif ()

if (DEFINED ENV{APHROS_PREFIX})
  set(CMAKE_INSTALL_PREFIX $ENV{APHROS_PREFIX}
      CACHE PATH "Install path prefix. Detected from environment APHROS_PREFIX" FORCE)
elseif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/prefix"
      CACHE PATH "Install path prefix. Defaults to 'prefix' in build directory." FORCE)
endif ()

set(BuildTypeValues None Debug Release RelWithDebInfo MinSizeRel)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: ${BuildTypeValues}." FORCE)
endif ()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BuildTypeValues})

set(CMAKEPREFIX ${CMAKE_INSTALL_PREFIX}/cmake/aphros)

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

if (DEFINED ENV{APHROS_PROFILE})
  set(Profile $ENV{APHROS_PROFILE}
      CACHE STRING "Profile for tools and wrappers, options are: ${ProfileValues}")
else()
  set(Profile "${ProfileDefault}"
      CACHE STRING "Profile for tools and wrappers, options are: ${ProfileValues}")
endif ()

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
    gitsnap
    hausdorff
    hist
    histbin
    intcol
    kedr
    kill
    line
    makearg
    makesim
    merge_hdf_xmf
    merge_hdf_xmf_all
    merge_raw_xmf
    merge_raw_xmf_all
    mpirun
    order
    part
    part0
    plain2vtk
    pngmp4
    py2conf
    run
    run0
    scale
    spawn
    split_hdf_xmf_all
    split_raw_xmf_all
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
      "${ScriptsDir}/plottools"
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
  if ("${name}" STREQUAL create_a_conf)
    set(deps ${deps} "${ScriptsDir}/sim_a.conf")
  endif()
  if ("${name}" STREQUAL create_base_conf)
    set(deps ${deps} "${ScriptsDir}/sim_base.conf")
  endif()
  set(output "${GenOutputDir}/${name}")
  list(APPEND GenOutputList "${output}")
  add_custom_command(
      OUTPUT "${output}"
      COMMAND ${CMAKE_COMMAND}
        -E env
        PREFIX=${CMAKE_INSTALL_PREFIX}
        PROFILE=${Profile}
        "${gen}" > "${output}"
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

if (USE_HYPRE)
  ExternalProject_Add(hypre
      URL "https://github.com/hypre-space/hypre/archive/v2.18.2.tar.gz"
      CMAKE_ARGS
        -DHYPRE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DHYPRE_SHARED=ON
        -DCMAKE_MACOSX_RPATH=OFF
        -DCMAKE_INSTALL_NAME_DIR=${CMAKE_INSTALL_PREFIX}/lib
      INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
      SOURCE_SUBDIR src
    )
endif()

if (USE_FPZIP)
    ExternalProject_Add(fpzip
        URL "https://github.com/LLNL/fpzip/releases/download/1.3.0/fpzip-1.3.0.tar.gz"
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_INSTALL_CMAKEDIR:STRING=cmake/fpzip
        -DBUILD_SHARED_LIBS=ON
        -DBUILD_TESTING=OFF
        -DBUILD_UTILITIES=OFF
        -DFPZIP_ENABLE_PIC=ON
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        )
endif ()

# AMGX
set(CUDA_ARCH 60 CACHE STRING "Cuda architecture version to be passed to AMGX")
set_property(CACHE CUDA_ARCH PROPERTY STRINGS 60 70)
if (USE_AMGX)
    ExternalProject_Add(amgx
        GIT_REPOSITORY "https://github.com/NVIDIA/AMGX"
        GIT_TAG v2.1.x
        GIT_SHALLOW ON
        PATCH_COMMAND git apply "${CMAKE_CURRENT_SOURCE_DIR}/amgx_no_comm_report.patch"
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DCUDA_ARCH=${CUDA_ARCH}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        )
endif ()
