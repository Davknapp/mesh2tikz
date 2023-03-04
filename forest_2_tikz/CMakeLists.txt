project(forest_2_tikz)

add_library(forest_2_tikz SHARED forest_to_tikz.cxx forest_to_tikz.hxx)

# Add t8code include directories
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE})
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/matrices)
target_link_libraries(forest_2_tikz transformation_matrix)
