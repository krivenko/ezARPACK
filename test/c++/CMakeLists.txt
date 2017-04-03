find_package(TriqsTest)

set(TESTS symmetric asymmetric complex)

foreach(t ${TESTS})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    set(t arpack_${t})
    add_executable(${t} ${s})
    target_link_libraries(${t} ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST} ${ARPACK_LIB})
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()
