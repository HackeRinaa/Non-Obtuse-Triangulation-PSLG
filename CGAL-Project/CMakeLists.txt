cmake_minimum_required(VERSION 3.1)  # Specify the minimum CMake version
project(TriangulationProject)        # Name of your project

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)           # Change to 14 or 17 if needed
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Enforce the standard

# Find and configure CGAL
find_package(CGAL REQUIRED)

# Include nlohmann/json library
include(FetchContent)
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3  # You can change this to the latest version
)
FetchContent_MakeAvailable(json)

# Specify the source files located in the /CodeFiles directory
set(SOURCE_FILES 
    CodeFiles/main.cpp 
    CodeFiles/triangulation_utils.cpp
)

# Add the executable that will be built from the source files
add_executable(triangulation_program ${SOURCE_FILES})

# Link the CGAL libraries and nlohmann/json to your executable
target_link_libraries(triangulation_program CGAL::CGAL nlohmann_json::nlohmann_json)

# Enable compiler warnings
if(MSVC)
    target_compile_options(triangulation_program PRIVATE /W4)  # For MSVC
else()
    target_compile_options(triangulation_program PRIVATE -Wall -Wextra -pedantic)  # For GCC/Clang
endif()

# Optionally, set installation rules
install(TARGETS triangulation_program DESTINATION bin)