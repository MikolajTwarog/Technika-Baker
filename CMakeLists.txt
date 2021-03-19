cmake_minimum_required(VERSION 3.17)
project(Technika_Baker)

set(CMAKE_CXX_STANDARD 14)

add_executable(Technika_Baker
        outerplanar/baker-outer-planar.hpp outerplanar/problems.hpp utils/level_face_traversal.h k-outerplanar/baker-k-outer-planar.hpp k-outerplanar/problems2.hpp teter.cpp)

add_subdirectory(Boost_tests)