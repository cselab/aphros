cmake_minimum_required(VERSION 3.3.0)

# Allow using <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

if (DEFINED ENV{APHROS_PREFIX})
  set(CMAKE_INSTALL_PREFIX $ENV{APHROS_PREFIX}
      CACHE PATH "Install path prefix. Detected from environment APHROS_PREFIX" FORCE)
else()
  message(FATAL_ERROR
    "Environment variable APHROS_PREFIX not set. Use `. ap.setenv`")
endif()

set(APHROS_PREFIX $ENV{APHROS_PREFIX})

# default build type
set(BuildTypeValues None Debug Release RelWithDebInfo MinSizeRel)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: ${BuildTypeValues}." FORCE)
endif ()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${BuildTypeValues})

# C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(T "openmp")
add_library(${T} INTERFACE IMPORTED)
if (APHROS_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  find_package(Threads REQUIRED)
  set_property(TARGET ${T} PROPERTY INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}")
  set_property(TARGET ${T} PROPERTY INTERFACE_LINK_LIBRARIES "${OpenMP_CXX_FLAGS}" Threads::Threads)
endif()

if (APHROS_USE_MPI)
  if ((NOT DEFINED MPI_C_COMPILER) AND (NOT DEFINED MPI_CXX_COMPILER))
    set(MPI_CXX_SKIP_MPICXX TRUE)
    find_package(MPI)
    if (NOT ${MPI_FOUND})
      message(FATAL_ERROR
        "**********\n"
        "MPI library wrapper is not found. Run cmake with -DUSE_MPI=0 to build aphros without MPI. Alternativly, you can install MPI by\n"
        "$ sudo apt install mpich\n"
        "**********\n")
     endif()
  endif()
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if (APHROS_USE_MPI AND APHROS_FIND_HDF)
  if (APHROS_USE_MPI)
    set(HDF5_PREFER_PARALLEL ON)
  endif()
  find_package(HDF5 COMPONENTS C HL)
  if (NOT ${HDF5_FOUND})
        message(FATAL_ERROR
        "**********\n"
        "HDF5 library is not found. Run cmake with \
        -DUSE_HDF=0 \
        -DUSE_BACKEND_CUBISM=0 \
        -DUSE_BACKEND_LOCAL=1 \
        -DUSE_BACKEND_NATIVE=1 \
        to build aphros without HDF5. Alternativly, you can install HDF5 by\n"
        "$ sudo apt install libhdf5-mpich-dev hdf5-tools\n"
        "**********\n")
  endif()
  if (APHROS_USE_MPI AND NOT HDF5_IS_PARALLEL)
        message(FATAL_ERROR
        "**********\n"
        "A serial HDF5 library is found, but aphros needs a parallel HDF5. Run cmake with \
        -DUSE_HDF=0 \
        -DUSE_BACKEND_CUBISM=0 \
        -DUSE_BACKEND_LOCAL=1 \
        -DUSE_BACKEND_NATIVE=1 \
        to build aphros without HDF5. Alternativly, you can install parallel HDF5 by\n"
        "$ sudo apt install libhdf5-mpich-dev hdf5-tools\n"
        "**********\n")
  endif()
endif()

# warnings
add_compile_options(-Wall -pedantic -Wextra)
