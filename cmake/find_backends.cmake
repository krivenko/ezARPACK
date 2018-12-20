message(STATUS "Detecting matrix algebra backends")

# Find Eigen3
find_package(Eigen3 3.3 CONFIG)
if(Eigen3_FOUND)
  message(STATUS "Found Eigen3 version ${Eigen3_VERSION}")
endif(Eigen3_FOUND)

# Find Blaze
find_package(blaze 3.0 QUIET CONFIG)
if(blaze_FOUND)
  message(STATUS "Found Blaze version ${blaze_VERSION}")
endif(blaze_FOUND)

# Armadillo
find_package(Armadillo CONFIG)

# Find TRIQS
find_package(Cpp2Py CONFIG)
find_package(TRIQS CONFIG)
