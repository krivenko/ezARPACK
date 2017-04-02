find_package(TriqsTest)

# Detect ARPACK-NG
find_package(ARPACK REQUIRED)

set(CMAKE_BUILD_TYPE Release)
add_definitions(${TRIQS_CXX_DEFINITIONS})
include_directories(${CMAKE_SOURCE_DIR}/c++ ${TRIQS_INCLUDE_ALL})

set(ARPACK_TESTS symmetric asymmetric complex)

foreach(t ${ARPACK_TESTS})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    set(t arpack_${t})
    add_executable(${t} ${s})
    target_link_libraries(${t} ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST} ${ARPACK_LIB})
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()
