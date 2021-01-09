cmake_minimum_required(VERSION 3.16)
project(project2)

set(CMAKE_CXX_STANDARD 14)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(project2 project2.cpp atoms.cpp atoms.h molecule.h utils.h utils.cpp)
target_link_libraries (project2 Eigen3::Eigen)
