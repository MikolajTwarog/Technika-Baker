//#include "outerplanar/test.hpp"
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <vector>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "k-outerplanar/baker-k-outer-planar.hpp"
#include "k-outerplanar/problems2.hpp"
#include <vector>

int main(int argc, char** argv) {
    Graph g6(5);
    add_edge(0,1,g6);
    add_edge(1,2,g6);
    add_edge(2,3,g6);
    add_edge(3, 0, g6);

    add_edge(4, 5, g6);
    add_edge(5, 6, g6);
    add_edge(6, 4, g6);

    add_edge(4, 0, g6);
    add_edge(4, 1, g6);
    add_edge(5, 2, g6);
    baker2<independent_set>(g6);
}

