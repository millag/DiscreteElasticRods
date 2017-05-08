cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(DiscreteElasticRods)

# Find Qt
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
# Find the QtWidgets library
find_package(Qt5Widgets)

# CML should be a submodule
set(CML_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/libs/cml-1_0_2")
# DLIB should be a submodule
set(DLIB_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/libs/dlib-18.9")

file(GLOB_RECURSE PROJECT_HEADERS "include/*.h" "include/*.hpp")
file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp")
set(PROJECT_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include")

include_directories(${CML_INCLUDE_DIRS})
include_directories(${DLIB_INCLUDE_DIRS})
include_directories(${PROJECT_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} "${PROJECT_HEADERS};${PROJECT_SOURCES}")

# Use the Widgets module from Qt 5.
target_link_libraries(${PROJECT_NAME} Qt5::Widgets)
