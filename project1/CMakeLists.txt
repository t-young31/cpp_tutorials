cmake_minimum_required(VERSION 3.16)
project(project1)

set(CMAKE_CXX_STANDARD 14)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(project1 project1.cpp atoms.cpp atoms.h)
target_link_libraries (project1 Eigen3::Eigen)
