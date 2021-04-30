cmake_minimum_required(VERSION 3.17)
project(Technika_Baker)

set(CMAKE_CXX_STANDARD 14)

add_executable(Technika_Baker
        outerplanar/baker-outer-planar.hpp outerplanar/problems.hpp utils/level_face_traversal.hpp k-outerplanar/baker-k-outer-planar.hpp k-outerplanar/problems2.hpp teter.cpp utils/cyclic_vector.hpp "Baker's Technique/bakers_technique.hpp" utils/random_graph.hpp utils/visitors.hpp utils/find_outer_face.hpp utils/name_levels.hpp)

add_subdirectory(Boost_tests)