cmake_minimum_required(VERSION 3.3.0)

cmake_minimum_required(VERSION 2.8)
project(overlap)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Header-only library.
include_directories(${PROJECT_SOURCE_DIR})

# Add compiler flags for C++11.
if(${CMAKE_VERSION} VERSION_LESS "3.1")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
        "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
    else()
        message(FATAL_ERROR
            "Unsupported compiler detected: ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Set up the Eigen3 library.
find_package(Eigen3 3.0.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(EXE "a")
add_executable(${EXE} a.cpp)
