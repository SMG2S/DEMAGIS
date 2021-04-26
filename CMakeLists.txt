cmake_minimum_required(VERSION 3.17)

project(Artificial-Matrices LANGUAGES C CXX)

find_package( Boost COMPONENTS serialization program_options REQUIRED)

add_subdirectory(wrappers)

add_executable(matgen.exe matgen.cpp)
target_link_libraries(matgen.exe PRIVATE blas_lapack_comm Boost::program_options )




