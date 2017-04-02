# Install C++ headers
install(FILES arpack.hpp
        arpack_worker.hpp
        worker_asymmetric.hpp
        worker_complex.hpp
        worker_symmetric.hpp
        DESTINATION include/triqs/arrays/arpack)

# Generate and install version.hpp
configure_file(version.hpp.in version.hpp)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.hpp DESTINATION include/triqs/arrays/arpack)
