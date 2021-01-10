cmake_minimum_required(VERSION 3.16)
project(project3)

set(CMAKE_CXX_STANDARD 14)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(project3 project3.cpp utils.h utils.cpp)
target_link_libraries (project3 Eigen3::Eigen)
