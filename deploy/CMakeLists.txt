cmake_minimum_required(VERSION 3.3.0)
project("aphros-deploy")

set(NAMEPREFIX ap.)

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

# Hypre
set(hypre_tar "${CMAKE_CURRENT_SOURCE_DIR}/extern/hypre.tar.gz")
set(hypre_url "https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.11.2.tar.gz")
ExternalProject_Add(hypre
    URL ${hypre_tar}
    CMAKE_ARGS
      -DHYPRE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DHYPRE_SHARED=ON
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    SOURCE_SUBDIR src
  )

# Eigen
set(eigen_tar "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen.tar.bz2")
set(eigen_url "http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2")
ExternalProject_Add(eigen
    URL ${eigen_tar}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:STRING=${CMAKE_INSTALL_PREFIX}
      -DINCLUDE_INSTALL_DIR:STRING=include
      -DCMAKEPACKAGE_INSTALL_DIR:STRING=cmake
      -DBUILD_TESTING=OFF
      -DEIGEN_BUILD_PKGCONFIG=OFF
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )

# fpzip
set(fpzip_tar "${CMAKE_CURRENT_SOURCE_DIR}/extern/fpzip.tar.gz")
set(fpzip_url "https://computing.llnl.gov/projects/floating-point-compression/download/fpzip-1.2.0.tar.gz")
ExternalProject_Add(fpzip
    URL ${fpzip_tar}
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_COMMAND make -C src CXXFLAGS=-O3\ -fPIC
    INSTALL_COMMAND ""
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )
set(fpzip_build "${CMAKE_CURRENT_BINARY_DIR}/fpzip-prefix/src/fpzip")
install(FILES ${fpzip_build}/lib/libfpzip.a DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${fpzip_build}/inc/fpzip.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

add_subdirectory(lib)
