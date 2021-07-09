cmake_minimum_required(VERSION 3.18)

project(Artificial-Matrices LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

find_package( Boost COMPONENTS serialization program_options REQUIRED)
find_package(MPI REQUIRED)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(wrappers)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/wrappers)

add_subdirectory(examples)

