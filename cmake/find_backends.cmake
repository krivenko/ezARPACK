message(STATUS "Detecting matrix algebra backends")

# Find Eigen3
find_package(Eigen3)

# Find TRIQS
find_package(Cpp2Py QUIET)
find_package(TRIQS QUIET)
