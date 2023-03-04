project(transformation_matrix)

# Add t8code include directories
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE})
INCLUDE_DIRECTORIES($ENV{T8_INCLUDE}/src)
add_library(transformation_matrix SHARED transformation_matrix.cxx transformation_matrix.hxx)


