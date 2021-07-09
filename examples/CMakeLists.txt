add_executable(driver_lapack.exe driver_lapack.cpp)
target_link_libraries(driver_lapack.exe PRIVATE blas_lapack_comm Boost::program_options )

add_executable(driver_scalapack.exe driver_scalapack.cpp)
target_link_libraries(driver_scalapack.exe PRIVATE blas_lapack_comm Boost::program_options MPI::MPI_CXX )

