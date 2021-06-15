cmake_minimum_required(VERSION 3.17)
project(Technika_Baker)

set(CMAKE_CXX_STANDARD 14)

add_executable(Technika_Baker
        utils/level_face_traversal.hpp baker/baker-k-outer-planar.hpp baker/problems.hpp utils/cyclic_vector.hpp "bakers_technique/bakers_technique.hpp" utils/random_graph.hpp utils/visitors.hpp utils/name_levels.hpp bodlaender/create_tree_decomposition.hpp bodlaender/bodlaender_impl.hpp example.cpp)

add_subdirectory(Boost_tests)