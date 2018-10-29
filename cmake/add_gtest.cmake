set(GOOGLETEST_VERSION "1.8.1")

function(add_gtest)
  if(NOT BUILD_GTEST)
    find_package(GTest ${GOOGLETEST_VERSION})
  endif(NOT BUILD_GTEST)

  if(BUILD_GTEST OR (NOT GTEST_FOUND))
    set(GTEST_TARBALL_URL "https://github.com/google/googletest/archive/release-${GOOGLETEST_VERSION}.tar.gz")
    set(GTEST_TARBALL "${CMAKE_BINARY_DIR}/googletest-${GOOGLETEST_VERSION}.tar.gz")

    if(NOT EXISTS "${GTEST_TARBALL}")
      message(STATUS "Downloading GTest")
      file(DOWNLOAD ${GTEST_TARBALL_URL}
           ${GTEST_TARBALL}
           STATUS GTEST_DOWNLOAD_STATUS
           SHOW_PROGRESS)
      list(GET GTEST_DOWNLOAD_STATUS 0 status_code)
      list(GET GTEST_DOWNLOAD_STATUS 1 status_msg)
      if(NOT status_code EQUAL 0)
        message(FATAL_ERROR "Could not download ${GTEST_TARBALL_URL}: ${status_msg}")
      endif(NOT status_code EQUAL 0)
    endif(NOT EXISTS "${GTEST_TARBALL}")

    if(NOT IS_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-release-${GOOGLETEST_VERSION}")
      execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${GTEST_TARBALL} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endif(NOT IS_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-release-${GOOGLETEST_VERSION}")

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    message(STATUS "Configuring GTest")
    add_subdirectory(${CMAKE_BINARY_DIR}/googletest-release-${GOOGLETEST_VERSION} "gtest" EXCLUDE_FROM_ALL)
  endif(BUILD_GTEST OR (NOT GTEST_FOUND))
endfunction()
