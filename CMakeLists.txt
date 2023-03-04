
#
# This file was copied and modified from https://cmake.org/examples/
#

# CMakeLists files in this project can
# refer to the root source directory of the project as ${PNG2MESH_SOURCE_DIR} and
# to the root binary directory of the project as ${PNG2MESH_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (PNG2MESH)
add_subdirectory(matrices)
add_subdirectory(forest_2_tikz)

# Set compilers to mpicc and mpicxx
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Add t8code include directories
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE})
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE}/src)
# Add t8code library directories
LINK_DIRECTORIES($ENV{T8_LIB})

#Add mesh2tikz library with source-files.
add_library (mesh2tikz mesh2tikz_demo.cxx forest_2_tikz/forest_to_tikz.hxx)

add_executable(mesh2tikz_demo mesh2tikz_demo.cxx forest_2_tikz/forest_to_tikz.hxx)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/matrices)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/forest_2_tikz)
#Link library against t8code
target_link_libraries(mesh2tikz t8 sc)

target_include_directories(mesh2tikz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/matrices)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/forest_2_tikz)
target_link_libraries(mesh2tikz_demo LINK_PUBLIC mesh2tikz transformation_matrix forest_2_tikz)


