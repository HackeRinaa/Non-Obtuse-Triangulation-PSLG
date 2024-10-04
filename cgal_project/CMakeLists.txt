cmake_minimum_required(VERSION 3.1)  # Specify the minimum CMake version
project(TriangulationProject)        # Name of your project

# Find and configure CGAL
find_package(CGAL REQUIRED)

# Add the executable that will be built from the source files
add_executable(triangulation_program main.cpp triangulation_utils.cpp)

# Link the CGAL libraries to your executable
target_link_libraries(triangulation_program CGAL::CGAL)
