cmake_minimum_required(VERSION 3.1)  # Specify the minimum CMake version
project(TriangulationProject)        # Name of your project

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)           # Change to 14 or 17 if needed
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Enforce the standard

# Specify the path to find Qt5
set(CMAKE_PREFIX_PATH "/usr/lib/x86_64-linux-gnu/cmake/Qt5")

# Find and configure CGAL
find_package(CGAL REQUIRED)

# Find and configure Qt
find_package(Qt5 COMPONENTS Widgets REQUIRED)

# Specify the source files located in the /CodeFiles directory
set(SOURCE_FILES 
    CodeFiles/main.cpp 
    CodeFiles/main_visualization.cpp 
    CodeFiles/triangulation_utils.cpp
    CodeFiles/visualization.cpp
)

# Add the executable that will be built from the source files
add_executable(triangulation_program ${SOURCE_FILES})

# Link the CGAL and Qt libraries to your executable
target_link_libraries(triangulation_program CGAL::CGAL Qt5::Widgets)

# Enable compiler warnings
if(MSVC)
    target_compile_options(triangulation_program PRIVATE /W4)  # For MSVC
else()
    target_compile_options(triangulation_program PRIVATE -Wall -Wextra -pedantic)  # For GCC/Clang
endif()

# Optionally, set installation rules
install(TARGETS triangulation_program DESTINATION bin)
