
#
# This file was copied and modified from https://cmake.org/examples/
#

# CMakeLists files in this project can
# refer to the root source directory of the project as ${PNG2MESH_SOURCE_DIR} and
# to the root binary directory of the project as ${PNG2MESH_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (PNG2MESH)

# Set compilers to mpicc and mpicxx
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Add t8code include directories
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE})
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE}/src)
# Add t8code library directories
LINK_DIRECTORIES($ENV{T8_LIB})


