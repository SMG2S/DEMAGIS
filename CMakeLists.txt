cmake_minimum_required(VERSION 3.18)

project(Artificial-Matrices LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

find_package( Boost COMPONENTS serialization program_options REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(wrappers)

add_executable(driver_lapack.exe driver_lapack.cpp)
target_link_libraries(driver_lapack.exe PRIVATE blas_lapack_comm Boost::program_options )

#add_executable(matgen_dist.exe matgen_dist.cpp)
#target_link_libraries(matgen_dist.exe PRIVATE blas_lapack_comm Boost::program_options MPI::MPI_CXX )



