# Start configuration
cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)

# Version number
project(triqs_arpack CXX)

set(TRIQS_ARPACK_VERSION "0.6")
message(STATUS "Configuring triqs_arpack version ${TRIQS_ARPACK_VERSION}")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${TRIQS_ROOT}/share/cmake)

# Load Cpp2Py and TRIQS
find_package(Cpp2Py 1.5 EXACT REQUIRED)
find_package(TRIQS 1.5 EXACT REQUIRED)

# Get hash
triqs_get_git_hash_of_source_dir(TRIQS_ARPACK_GIT_HASH)
message(STATUS "Git hash: ${TRIQS_ARPACK_GIT_HASH}")

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT})

# Install FindARPACK.cmake
install(FILES cmake/FindARPACK.cmake DESTINATION share/cmake)

option(Tests "Enable Tests" ON)
option(Examples "Build examples" ON)

add_compile_options(-std=c++14)

# Do we need to build anything?
if(Tests OR Examples)
    # Detect ARPACK-NG
    find_package(ARPACK REQUIRED)

    set(CMAKE_BUILD_TYPE Release)
    add_definitions(${TRIQS_CXX_DEFINITIONS})
    include_directories(${CMAKE_SOURCE_DIR}/c++ ${TRIQS_INCLUDE_ALL})
endif(Tests OR Examples)

# Install C++ headers
add_subdirectory(c++)

# Build automated tests
if(Tests)
    enable_testing()
    add_subdirectory(test)
endif(Tests)

# Build examples
if(Examples)
    add_subdirectory(example)
endif(Examples)
