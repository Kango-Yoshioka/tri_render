cmake_minimum_required(VERSION 3.22)
project(tri_render)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES src/main.cpp src/Camera.cpp src/Camera.h src/TriMesh.cpp src/TriMesh.h src/util.cpp src/util.h)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

include("C:/CGProgramming/Libraries/use_opengl_eigen.cmake")