message(STATUS "Detecting matrix algebra backends")

# Find Eigen3
find_package(Eigen3 CONFIG)
if(Eigen3_FOUND)
  if(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${Eigen3_VERSION}")
    if(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
      add_compile_definitions(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
  else(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${EIGEN3_VERSION_STRING}")
    if(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
      add_compile_definitions(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
  endif(Eigen3_VERSION)
endif(Eigen3_FOUND)

# Find Blaze
find_package(blaze 3.0 QUIET CONFIG)
if(blaze_FOUND)
  message(STATUS "Found Blaze version ${blaze_VERSION}")
endif(blaze_FOUND)

# Armadillo
find_package(Armadillo CONFIG)

# Boost uBLAS
find_package(Boost 1.58)

# Find TRIQS
find_package(Cpp2Py CONFIG)
find_package(TRIQS CONFIG)
