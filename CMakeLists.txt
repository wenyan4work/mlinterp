project(MLINTERP)
cmake_minimum_required(VERSION 3.10)

# C++11
set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP REQUIRED)

# mlinterp
include_directories(${PROJECT_SOURCE_DIR}/mlinterp)

# Set appropriate optimization flags
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(examples)
