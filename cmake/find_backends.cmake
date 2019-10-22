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
find_package(Armadillo QUIET CONFIG)
if(NOT Armadillo_FOUND)
  if(Armadillo_ROOT)
    find_library(ARMADILLO_LIBRARY NAMES armadillo PATHS ${Armadillo_ROOT}/lib NO_DEFAULT_PATH)
  else(Armadillo_ROOT)
    find_library(ARMADILLO_LIBRARY NAMES armadillo)
  endif(Armadillo_ROOT)

  if(ARMADILLO_LIBRARY)
    set(Armadillo_FOUND TRUE)
    get_filename_component(ARMADILLO_PREFIX ${ARMADILLO_LIBRARY} DIRECTORY)
    add_library(armadillo INTERFACE)
    target_include_directories(armadillo INTERFACE ${ARMADILLO_PREFIX}/include)
    target_link_libraries(armadillo INTERFACE ${ARMADILLO_LIBRARY})
  endif(ARMADILLO_LIBRARY)
endif(NOT Armadillo_FOUND)
if(Armadillo_FOUND)
  message(STATUS "Found Armadillo")
endif(Armadillo_FOUND)

# Boost uBLAS
find_package(Boost 1.58)

# Find TRIQS
find_package(Cpp2Py CONFIG)
find_package(TRIQS CONFIG)

# Find xtensor
find_package(xtensor CONFIG 0.20)
find_package(xtensor-blas CONFIG 0.16)
if(xtensor_FOUND AND xtensor-blas_FOUND)
  message(STATUS "Found xtensor version ${xtensor_VERSION}")
  message(STATUS "Found xtensor-blas version ${xtensor-blas_VERSION}")
endif(xtensor_FOUND AND xtensor-blas_FOUND)
