cmake_minimum_required(VERSION 3.10)
project(QuantumPIMD)

set(CMAKE_CXX_STANDARD 17)

# 1. Find Dependencies
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

# 2. Include paths
include_directories(
    ${PROJECT_SOURCE_DIR}/include 
    ${PROJECT_SOURCE_DIR}/external/glad/include
)

# 3. Explicit Source List
# We must include the .c file for GLAD here!
set(SOURCES 
    src/main.cpp
    src/engine/PIMDSystem.cpp
    external/glad/src/glad.c
)

add_executable(QuantumPIMD ${SOURCES})

# 4. Linking
# On Windows, we often need to link opengl32 explicitly for glad
target_link_libraries(QuantumPIMD glfw glm::glm opengl32)