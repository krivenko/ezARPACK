set(EXAMPLES example)

foreach(t ${EXAMPLES})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    add_executable(${t} ${s})
    target_link_libraries(${t} ${TRIQS_LIBRARY_ALL} ${ARPACK_LIB})
    triqs_set_rpath_for_target(${t})
endforeach()
